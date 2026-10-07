""" Web viewer for sequencing dots and cells of a finished pipeline run.

    python -m viewer.server /path/to/run --port 8000

then tunnel the port to your computer (see viewer/README.md) and open http://localhost:8000
"""

import os
import io
import re
import sys
import json
import time
import secrets
import argparse
import threading
from typing import List

import numpy as np
from PIL import Image
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import FileResponse, JSONResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles
from starlette.middleware.trustedhost import TrustedHostMiddleware

from .imaging import RunImages, load_config, find_input_files
from .tables import RunTables, parse_tile
from .files import LAYOUTS, MissingFiles, RunFiles, well_names


STATIC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static')
MAX_CROP = 2048
MAX_LAYERS = 128
TOKEN_FILE = os.path.expanduser('~/.config/starcall-viewer/token')
# the page and its scripts only load from the server itself (plus the Google font); images may be
# data: or blob: URLs, which the PNG/SVG export uses
CONTENT_SECURITY_POLICY = ("default-src 'self'; img-src 'self' data: blob:; style-src 'self' https://fonts.googleapis.com; "
                           "font-src https://fonts.gstatic.com; object-src 'none'; base-uri 'none'; frame-ancestors 'none'")


class Run:
    """ A pipeline run: its wells, and a Viewer for each well once its files are known to exist """

    def __init__(self, run_dir, cache_dir=None, layout='auto', **table_names):
        self.run_dir = os.path.abspath(run_dir)
        self.config = load_config(self.run_dir)
        self.inputfiles = find_input_files(self.run_dir, self.config)
        self.wells = well_names(self.config, self.inputfiles)
        # layout: 'tiled', 'untiled', or 'auto' to detect it for each well from the tables present
        self.files = RunFiles(self.run_dir, self.config, self.inputfiles, layout=layout, **table_names)
        self.cache_dir = cache_dir
        self.table_names = table_names
        self.viewers = {}
        self.lock = threading.Lock()

    def viewer(self, well):
        if well not in self.wells:
            raise HTTPException(404, 'Unknown well {}; this run has {}'.format(well, ', '.join(self.wells)))
        with self.lock:
            if well not in self.viewers:
                self.files.require(self.files.well_files(well), 'Files needed to show images of {} are missing'.format(well))
                self.viewers[well] = Viewer(self.run_dir, well, self.config, self.inputfiles[well], self.files, self.files.layout(well),
                        cache_dir=self.cache_dir, **self.table_names)
            return self.viewers[well]

    def default_well(self):
        """ The first well whose files are all present, otherwise the first well """
        for well in self.wells:
            if all(f['exists'] for f in self.files.well_files(well)):
                return well
        return self.wells[0] if self.wells else None


def phenotyping_channel_lists(config):
    """ config phenotyping_channels as one list per phenotype cycle, as workflow/rules/config.smk does """
    channels = config['phenotyping_channels']
    if not isinstance(channels[0], list):
        channels = [channels]
    return [list(c) for c in channels]


def all_phenotyping_channels(config):
    """ Every phenotype channel name, in config order, without repeats """
    return list(dict.fromkeys(c for channels in phenotyping_channel_lists(config) for c in channels))


def segmentation_channels(config):
    """ segmentation: channels as names; an index refers to the first phenotype cycle's channels """
    first = phenotyping_channel_lists(config)[0]
    channels = config.get('segmentation', {}).get('channels', first[:2])
    return [c if isinstance(c, str) else first[c] for c in channels]


class Viewer:
    def __init__(self, run_dir, well, config, paths, files, layout, cache_dir=None, **table_names):
        self.run_dir = run_dir
        self.well = well
        self.config = config
        self.files = files
        self.scale = self.config['phenotype_scale'] // self.config['bases_scale']
        self.images = RunImages(self.run_dir, well, self.config, paths=paths)
        self.tables = RunTables(self.run_dir, well, self.config, layout, cache_dir=cache_dir, **table_names)

        self.sequencing_channels = list(self.config['sequencing_channels'])
        self.phenotyping_channel_lists = phenotyping_channel_lists(self.config)
        self.phenotyping_channels = all_phenotyping_channels(self.config)
        self.segmentation_channels = segmentation_channels(self.config)

        self._contrast = None
        self._contrast_lock = threading.Lock()
        self._masks = {}

    ##### contrast #####

    def contrast(self):
        """ Default display range per channel, shared by all cycles so they can be compared.
        1st and 99.9th percentiles of pixels sampled from a few frames of a few cycles.
        """
        with self._contrast_lock:
            if self._contrast is not None:
                return self._contrast
            path = os.path.join(self.tables.cache_dir, '{}_contrast.json'.format(self.well))
            if os.path.exists(path):
                with open(path) as ifile:
                    self._contrast = json.load(ifile)
                return self._contrast

            result = {}
            for kind, cycles in (('seq', self.images.seq_cycles), ('pt', self.images.pt_cycles)):
                if not cycles: continue
                # phenotype cycles can each image different channels, so sample all of them
                picks = cycles if kind == 'pt' else sorted({cycles[0], cycles[len(cycles) // 2], cycles[-1]})
                samples = {}
                for cycle in picks:
                    pixels = self.images.cycles[cycle].sample_pixels()
                    for i, name in enumerate(self.cycle_channels(cycle)[:pixels.shape[0]]):
                        samples.setdefault(name, []).append(pixels[i])
                result[kind] = {}
                for name, parts in samples.items():
                    values = np.concatenate(parts)
                    values = values[np.isfinite(values)]
                    if not values.size: continue
                    lo, hi = np.percentile(values, [1, 99.9])
                    result[kind][name] = {'vmin': float(lo), 'vmax': float(max(hi, lo + 1)), 'max': float(values.max())}
            with open(path, 'w') as ofile:
                json.dump(result, ofile, indent=1)
            self._contrast = result
            return result

    def cycle_channels(self, cycle):
        """ Channel names of a cycle. Phenotype cycles PT, P1, P2... use the matching list of
        phenotyping_channels, or the last one if there are fewer lists than cycles.
        """
        if not self.images.is_pt(cycle):
            return self.sequencing_channels
        lists = self.phenotyping_channel_lists
        match = re.fullmatch(r'P(\d+)', cycle)
        index = 0 if cycle == 'PT' else int(match.group(1)) if match else self.images.pt_cycles.index(cycle)
        return lists[min(index, len(lists) - 1)]

    def channel_index(self, cycle, channel):
        channels = self.cycle_channels(cycle)
        if channel not in channels:
            raise HTTPException(400, 'Unknown channel {} for cycle {}'.format(channel, cycle))
        return channels.index(channel)

    ##### views #####

    def view_box(self, center=None, cell_box_pt=None, size=64, pad=0.25, min_size=32):
        """ A square view box in bases pixels of the well composite, either `size` pixels around
        center (bases pixels), or containing a cell bbox (phenotype pixels) with some padding.
        """
        if cell_box_pt is not None:
            r0, c0, r1, c1 = np.asarray(cell_box_pt, float) / self.scale
            side = max(r1 - r0, c1 - c0) * (1 + 2 * pad)
            side = int(np.ceil(max(side, min_size)))
            center = ((r0 + r1) / 2, (c0 + c1) / 2)
        else:
            side = int(size)
        r0 = int(np.floor(center[0] - side / 2))
        c0 = int(np.floor(center[1] - side / 2))
        return [r0, c0, r0 + side, c0 + side]

    def require_tile(self, tile, kinds):
        self.files.require(self.files.tile_files(self.well, tile, kinds),
                'Tables needed to show {} {} are missing'.format(self.well, tile))

    def check_tile(self, tile):
        """ The canonical name of a tile of this well's grid ('tile02x03y', or 'whole' for an untiled
        run); 404 for any other """
        if tile not in self.tables.sections:
            try:
                tile, _, _ = parse_tile(tile)
            except ValueError as error:
                raise HTTPException(404, str(error))
        if tile not in self.tables.sections:
            raise HTTPException(404, 'No tile ' + tile)
        return tile

    def dot(self, tile, index, size=None):
        """ size: None fits the view to the dot's cell (64 px around a dot with no cell),
        otherwise a size-pixel square centred on the dot
        """
        tile = self.check_tile(tile)
        self.require_tile(tile, ['dots'])
        rows = self.tables.dots(tile, [index])
        if len(rows) == 0:
            raise HTTPException(404, 'Dot {} not found in {} {}'.format(index, tile, self.tables.dots_table))
        row = rows.iloc[0]
        section = self.tables.section(tile)
        position = [float(row.position_x) + section[0], float(row.position_y) + section[1]]

        cell = None
        label = int(row.cell)
        if label > 0:
            self.require_tile(tile, ['dots', 'cells'])
            cell_id, cell_row = self.tables.cell_by_label(tile, label)
            cell_box = self.tables.cell_box(tile, cell_row)
            box = self.view_box(cell_box_pt=cell_box) if size is None else self.view_box(center=position, size=size)
            cell = {'id': cell_id, 'label': label, 'bbox_pt': [int(v) for v in cell_box]}
            cell['outline'] = self.outline(tile, label, cell_row, [int(v) * self.scale for v in box])
            cell['nucleus_outline'] = self.nucleus_outline(tile, label, [int(v) * self.scale for v in box])
        else:
            box = self.view_box(center=position, size=64 if size is None else size)

        cycles = self.tables.num_cycles(rows)
        bases = list(dict.fromkeys(c.split('_')[-1] for c in rows.columns if c.startswith('values_cycle')))
        return {
            'kind': 'dot',
            'tile': tile,
            'index': int(index),
            'sequence': str(row.max_seq),
            'cycles': cycles,
            'chastity': [float(row['chastity_cycle' + c]) for c in cycles],
            'mean_chastity': float(row.get('mean_chastity', np.nan)),
            'min_chastity': float(row.get('min_chastity', np.nan)),
            'zscore_bases': bases,
            'zscores': {b: [float(row['values_cycle{}_{}'.format(c, b)]) for c in cycles] for b in bases},
            'position': position,
            'cell': cell,
            **self.box_info(tile, box),
            'dots': [self.marker(index, position, box, row, current=True)],
        }

    def cell(self, cell_id, size=None):
        """ size: None fits the view to the cell, otherwise a size-pixel square centred on the cell """
        try:
            tile, label, row = self.tables.cell_by_id(cell_id)
        except KeyError as error:
            # the cell may be in a tile whose cell table is missing
            files = [f for tile in self.tables.tiles() for f in self.files.tile_files(self.well, tile, ['cells'])]
            if not all(f['exists'] for f in files):
                raise MissingFiles('Cell {} was not found in {}, and some cell tables are missing'.format(cell_id, self.well), files)
            raise HTTPException(404, str(error))
        self.require_tile(tile, ['cells', 'reads', 'dots'])
        cell_box = self.tables.cell_box(tile, row)
        if size is None:
            box = self.view_box(cell_box_pt=cell_box)
        else:
            center = ((cell_box[0] + cell_box[2]) / 2 / self.scale, (cell_box[1] + cell_box[3]) / 2 / self.scale)
            box = self.view_box(center=center, size=size)
        reads = self.tables.cell_reads(tile, cell_id)

        if 'dot_indicies' in reads:
            dot_source = 'dot_indicies'
            dots = self.tables.dots(tile, reads['dot_indicies'].split(':'))
        else:
            # cells_reads.csv made before the dot_indicies column existed
            dot_source = 'cell column of {} (no dot_indicies in {})'.format(self.tables.dots_table, self.tables.reads_table)
            dots = self.tables.dots_in_cell(tile, label)
        section = self.tables.section(tile)

        markers = []
        for index, dot in dots.iterrows():
            position = [float(dot.position_x) + section[0], float(dot.position_y) + section[1]]
            markers.append(self.marker(index, position, box, dot))

        read_list = []
        i = 0
        while 'read_{}'.format(i) in reads:
            read_list.append({key: reads.get('{}_{}'.format(key, i)) for key in
                    ('read', 'count', 'chastities', 'barcode_matches', 'barcode_hamming_dist')})
            i += 1

        return {
            'kind': 'cell',
            'tile': tile,
            'cell': {'id': int(cell_id), 'label': int(label), 'bbox_pt': [int(v) for v in cell_box],
                     'outline': self.outline(tile, label, row, [int(v) * self.scale for v in box]),
                     'nucleus_outline': self.nucleus_outline(tile, label, [int(v) * self.scale for v in box])},
            'num_reads': reads.get('num_reads'),
            'total_count': reads.get('total_count'),
            'reads': read_list,
            'dot_source': dot_source,
            **self.box_info(tile, box),
            'dots': markers,
        }

    def box_info(self, tile, box):
        box_pt = [int(v) * self.scale for v in box]
        return {
            'box': [int(v) for v in box],
            'box_pt': box_pt,
            'section': [int(v) for v in self.tables.section(tile)],
            # fraction of the crop each cycle imaged; 0 means no image data there in that cycle
            'coverage': {cycle: round(images.coverage(box_pt if self.images.is_pt(cycle) else box), 4)
                         for cycle, images in self.images.cycles.items()},
        }

    def marker(self, index, position, box, row, current=False):
        return {
            'index': int(index),
            'rel': [(position[0] - box[0]) / (box[2] - box[0]), (position[1] - box[1]) / (box[3] - box[1])],
            'sequence': str(row.max_seq),
            'min_chastity': float(row.min_chastity) if 'min_chastity' in row else None,
            'current': current,
        }

    ##### images #####

    def scaled_crop(self, tile, cycle, channel, box, vmin=None, vmax=None, contrast='global'):
        """ One channel of a crop scaled from vmin (0) to vmax (1), and a mask of the pixels some
        frame of the cycle imaged. box is in the cycle's own pixels.
        """
        if cycle not in self.images.cycles:
            raise HTTPException(400, 'Unknown cycle ' + cycle)
        pt = self.images.is_pt(cycle)
        index = self.channel_index(cycle, channel)
        image = self.images.crop(cycle, box, self.tables.section(tile, phenotype=pt))[index]

        if vmin is None or vmax is None:
            if contrast == 'local' and np.isfinite(image).any():
                lo, hi = np.nanpercentile(image, [0.5, 99.8])
            else:
                default = self.contrast()['pt' if pt else 'seq'][channel]
                lo, hi = default['vmin'], default['vmax']
            vmin = lo if vmin is None else vmin
            vmax = hi if vmax is None else vmax
        imaged = np.isfinite(image)
        scaled = np.clip((image - vmin) / max(vmax - vmin, 1e-6), 0, 1)
        return np.nan_to_num(scaled, nan=0), imaged

    def crop_png(self, tile, cycle, channel, box, vmin=None, vmax=None, contrast='global', color=None):
        """ 8-bit RGBA PNG of one channel of a crop, displayed from vmin (black) to vmax (full colour).
        color is a hex RGB string for a black -> colour map; without it the image is grayscale.
        Pixels outside every frame of the cycle (not imaged) are fully transparent.
        """
        tile = self.check_tile(tile)
        scaled, imaged = self.scaled_crop(tile, cycle, channel, parse_box(box), vmin, vmax, contrast)
        # pixels no frame imaged in this cycle (NaN) are transparent, so the page can mark them
        rgba = np.concatenate([scaled[..., None] * hex_rgb(color), imaged[..., None] * 255.0], axis=2)
        return png_bytes(Image.fromarray(rgba.astype(np.uint8), mode='RGBA'))

    def composite_png(self, tile, box, layers, contrast='global'):
        """ 8-bit RGBA PNG of several channels of a crop added together, each from black to its colour,
        like the channels of a multichannel tif (make_variant_cell_images) shown as a composite.
        box is in bases pixels. Each layer is 'cycle|channel|rrggbb' or 'cycle|channel|rrggbb|vmin|vmax'.
        The image is at phenotype resolution if any layer is a phenotype cycle (sequencing layers are
        enlarged to match), otherwise at sequencing resolution. Pixels no layer imaged are transparent.
        """
        tile = self.check_tile(tile)
        box = parse_box(box)
        layers = list(dict.fromkeys(layers))   # a layer ticked twice would just be twice as bright
        if not layers:
            raise HTTPException(400, 'No layers')
        if len(layers) > MAX_LAYERS:
            raise HTTPException(400, 'At most {} layers'.format(MAX_LAYERS))
        parsed = []
        for layer in layers:
            parts = layer.split('|')
            if len(parts) not in (3, 5) or not re.fullmatch('[0-9a-fA-F]{6}', parts[2]):
                raise HTTPException(400, 'Bad layer ' + layer)
            try:
                vmin, vmax = (float(parts[3]), float(parts[4])) if len(parts) == 5 else (None, None)
            except ValueError:
                raise HTTPException(400, 'Bad layer ' + layer)
            parsed.append((parts[0], parts[1], parts[2], vmin, vmax))

        scale = self.scale if any(self.images.is_pt(cycle) for cycle, *_ in parsed) else 1
        shape = ((box[2] - box[0]) * scale, (box[3] - box[1]) * scale)
        total = np.zeros(shape + (3,), np.float32)
        imaged = np.zeros(shape, bool)
        for cycle, channel, color, vmin, vmax in parsed:
            pt = self.images.is_pt(cycle)
            layer_box = [v * self.scale for v in box] if pt else box
            scaled, layer_imaged = self.scaled_crop(tile, cycle, channel, parse_box(layer_box), vmin, vmax, contrast)
            if not pt and scale > 1:
                # box_pt is box * scale, so each sequencing pixel covers exactly scale x scale phenotype pixels
                scaled = np.repeat(np.repeat(scaled, scale, axis=0), scale, axis=1)
                layer_imaged = np.repeat(np.repeat(layer_imaged, scale, axis=0), scale, axis=1)
            total += scaled[..., None] * hex_rgb(color)
            imaged |= layer_imaged
        rgba = np.concatenate([np.clip(total, 0, 255), imaged[..., None] * 255.0], axis=2)
        return png_bytes(Image.fromarray(rgba.astype(np.uint8), mode='RGBA'))

    def crop_labels(self, image, tile, box_pt):
        """ Crop of a tile-local segmentation label image (phenotype pixels) over box_pt, 0 outside it """
        section = self.tables.section(tile, phenotype=True)
        shape = (box_pt[2] - box_pt[0], box_pt[3] - box_pt[1])
        labels = np.zeros(shape, image.dtype)
        r0, c0 = box_pt[0] - section[0], box_pt[1] - section[1]
        q0, s0 = max(r0, 0), max(c0, 0)
        q1, s1 = min(r0 + shape[0], image.shape[0]), min(c0 + shape[1], image.shape[1])
        if q0 < q1 and s0 < s1:
            labels[q0 - r0:q1 - r0, s0 - c0:s1 - c0] = image[q0:q1, s0:s1]
        return labels

    def object_mask_in_box(self, kind, tile, label, row, box_pt):
        """ Boolean mask over box_pt (phenotype pixels of the well composite) of a cell (kind='cells')
        or of its matched nucleus (kind='nuclei', row from nuclei.csv)
        """
        full_mask = self.full_mask(tile, kind)
        if full_mask is not None:
            # exact segmentation, while the segmentation mask image still exists; the unmatched
            # nuclei mask is labelled by orig_index rather than by the matched cell label
            return self.crop_labels(full_mask, tile, box_pt) == (label if kind == 'cells' else row.orig_index)
        shape = (box_pt[2] - box_pt[0], box_pt[3] - box_pt[1])
        mask = np.zeros(shape, bool)
        object_mask = self.tables.cell_mask(row)
        object_box = self.tables.cell_box(tile, row)
        if object_mask is None:
            object_mask = np.ones((object_box[2] - object_box[0], object_box[3] - object_box[1]), bool)
        r0, c0 = object_box[0] - box_pt[0], object_box[1] - box_pt[1]
        q0, s0 = max(r0, 0), max(c0, 0)
        q1, s1 = min(r0 + object_mask.shape[0], shape[0]), min(c0 + object_mask.shape[1], shape[1])
        if q0 < q1 and s0 < s1:
            mask[q0:q1, s0:s1] = object_mask[q0 - r0:q1 - r0, s0 - c0:s1 - c0]
        return mask

    def outline(self, tile, label, row, box_pt, kind='cells'):
        """ Outline of a cell or nucleus over box_pt as an SVG path along the mask's pixel edges, in
        phenotype pixels of the box (x = column, y = row), so the browser can draw it as a thin line at any zoom.
        """
        return edge_path(self.object_mask_in_box(kind, tile, label, row, box_pt))

    def nucleus_outline(self, tile, label, box_pt):
        """ Outline of the nucleus matched to a cell, or None if there is no nuclei table """
        row = self.tables.nucleus_by_label(tile, label)
        return None if row is None else self.outline(tile, label, row, box_pt, kind='nuclei')

    def labels_in_box(self, tile, box_pt, exclude=None, kind='cells'):
        """ Tile-local cell labels over box_pt (0 = none) of cells or of their matched nuclei,
        leaving out the label exclude
        """
        full_mask = self.full_mask(tile, kind)
        if full_mask is not None:
            labels = self.crop_labels(full_mask, tile, box_pt).astype(np.int64)
            if kind == 'nuclei':
                # unmatched nuclei mask labels (orig_index) -> matched cell labels; unmatched nuclei become 0
                orig = self.tables.nucleus_table(tile).orig_index.to_numpy()
                lookup = np.zeros(max(int(orig.max()), int(labels.max())) + 1, np.int64)
                lookup[orig] = np.arange(1, len(orig) + 1)
                labels = lookup[labels]
            labels = labels.astype(np.int32)
        else:
            # without the mask image, paint each downscaled mask from the table (finest available)
            labels = np.zeros((box_pt[2] - box_pt[0], box_pt[3] - box_pt[1]), np.int32)
            labels_found, rows = self.tables.cells_in_box(tile, box_pt, kind)
            for label, (_, row) in zip(labels_found, rows.iterrows()):
                labels[self.object_mask_in_box(kind, tile, label, row, box_pt)] = label
        if exclude:
            labels[labels == exclude] = 0
        return labels

    def region(self, tile, box, exclude_cell=None):
        """ Every cell outline, nucleus outline and dot inside a view box (bases pixels), for the
        'all in view' option. Only this tile's tables are used, so dots and cells of a neighbouring
        tile are not included.
        """
        tile = self.check_tile(tile)
        self.require_tile(tile, ['dots', 'cells'])
        box = parse_box(box)
        box_pt = [v * self.scale for v in box]
        labels = self.labels_in_box(tile, box_pt, exclude=exclude_cell)
        nuclei = None
        if self.tables.has_nuclei(tile):
            nuclei = edge_path(self.labels_in_box(tile, box_pt, exclude=exclude_cell, kind='nuclei'))
        section = self.tables.section(tile)
        dots = self.tables.dots_in_box(tile, box)
        markers = []
        for index, dot in dots.iterrows():
            position = [float(dot.position_x) + section[0], float(dot.position_y) + section[1]]
            markers.append(self.marker(index, position, box, dot))
        # each cell separately too, so the enlarged view can open a cell when its outline is clicked
        ids = self.tables._cell_table(tile).index
        cells = [{'id': int(ids[label - 1]), 'label': int(label), 'outline': edge_path(labels == label),
                  'area': area_path(labels == label)}
                 for label in np.unique(labels[labels != 0])]
        return {
            'box': box,
            'outline': edge_path(labels),
            'nuclei_outline': nuclei,
            'num_cells': len(cells),
            'cells': cells,
            'dots': markers,
        }

    MASK_IMAGES = {
        # matched cell labels
        'cells': ['cells_mask.tif'],
        # labelled by nuclei.csv orig_index
        'nuclei': ['nuclei_mask_unmatched_grid{grid}.tif', 'nuclei_mask_unmatched.tif'],
    }

    def full_mask(self, tile, kind='cells'):
        """ Memory mapped segmentation label image of a tile, or None once it has been deleted """
        key = (tile, kind)
        if key not in self._masks:
            mask = None
            for name in self.MASK_IMAGES[kind]:
                path = self.tables.segmentation_path(tile, name.format(grid=self.tables.grid_size))
                if os.path.exists(path):
                    try:
                        import tifffile
                        mask = tifffile.memmap(path, mode='r')
                        break
                    except Exception:
                        mask = None
            self._masks[key] = mask
        return self._masks[key]

    def warm_up(self):
        """ Open nd2 files and compute default contrast in the background """
        def run():
            try:
                for cycle in self.images.cycles.values():
                    cycle.file() if cycle.path.endswith('.nd2') else None
                self.contrast()
            except Exception as error:
                print('warm up failed:', error, file=sys.stderr)
        threading.Thread(target=run, daemon=True).start()


def edge_path(labels):
    """ SVG path along every pixel edge between different values of a mask or label image,
    in pixels of the image (x = column, y = row)
    """
    # replicate the border so a cell cut off by the crop isn't outlined along the crop edge
    padded = np.pad(labels, 1, mode='edge')
    horizontal = padded[1:, 1:-1] != padded[:-1, 1:-1]   # edge above row y, shape (rows + 1, cols)
    vertical = padded[1:-1, 1:] != padded[1:-1, :-1]     # edge left of column x, shape (rows, cols + 1)
    parts = []
    for y, line in enumerate(horizontal):
        for start, stop in runs(line):
            parts.append('M{} {}H{}'.format(start, y, stop))
    for x, line in enumerate(vertical.T):
        for start, stop in runs(line):
            parts.append('M{} {}V{}'.format(x, start, stop))
    return {'size': [int(labels.shape[0]), int(labels.shape[1])], 'path': ''.join(parts)}


def area_path(mask):
    """ SVG path filling the pixels of a mask (one rectangle per run of pixels in each row), so the
    browser can tell when a click lands inside the cell
    """
    parts = []
    for y in np.nonzero(mask.any(axis=1))[0]:
        for start, stop in runs(mask[y]):
            parts.append('M{} {}h{}v1h-{}z'.format(start, y, stop - start, stop - start))
    return ''.join(parts)


def runs(line):
    """ (start, stop) of each run of True values in a 1d boolean array """
    edges = np.diff(np.concatenate([[0], line.astype(np.int8), [0]]))
    return zip(np.nonzero(edges == 1)[0], np.nonzero(edges == -1)[0])


def parse_box(box):
    """ A crop box (r0, c0, r1, c1), from a list or from 'r0,c0,r1,c1' """
    try:
        box = [int(v) for v in (box.split(',') if isinstance(box, str) else box)]
    except ValueError:
        raise HTTPException(400, 'Bad box')
    if len(box) != 4 or box[2] <= box[0] or box[3] <= box[1] or max(box[2] - box[0], box[3] - box[1]) > MAX_CROP:
        raise HTTPException(400, 'Bad box')
    return box


def hex_rgb(color):
    """ 'rrggbb' -> float RGB 0-255; white without a colour """
    return np.array([int(color[i:i + 2], 16) for i in (0, 2, 4)] if color else [255, 255, 255], np.float32)


def png_bytes(image):
    buffer = io.BytesIO()
    image.save(buffer, format='PNG', compress_level=1)
    return buffer.getvalue()


def default_hosts():
    """ Host names the browser may reach the server by: localhost through a tunnel, or this node """
    import socket
    return ['localhost', '127.0.0.1', socket.gethostname(), socket.getfqdn()]


def make_app(run, token, hosts=None):
    """ token is required on every request, as ?token= once (then a cookie). hosts are the Host
    header values accepted, which stops DNS rebinding (default: localhost and this node).
    """
    if not token:
        raise ValueError('The viewer needs an access token')
    app = FastAPI(title='starcall viewer', docs_url=None, redoc_url=None, openapi_url=None)

    @app.exception_handler(MissingFiles)
    def missing_files(request: Request, error: MissingFiles):
        return JSONResponse(error.as_dict(), status_code=424)

    # middleware added later runs first: host check, then headers, then the token
    @app.middleware('http')
    async def check_token(request: Request, call_next):
        from_query = request.query_params.get('token')
        given = from_query or request.cookies.get('viewer_token')
        if not secrets.compare_digest((given or '').encode(), token.encode()):
            return Response('Forbidden: open the link printed by the server (with ?token=...)', status_code=403)
        if from_query and request.url.path == '/':
            # the cookie is enough from here on; take the token out of the address bar and history
            # (the browser keeps the #dot/... or #cell/... part across the redirect)
            response = RedirectResponse('/', status_code=303)
        else:
            response = await call_next(request)
        response.set_cookie('viewer_token', token, max_age=30 * 24 * 3600, httponly=True, samesite='strict')
        return response

    @app.middleware('http')
    async def headers(request: Request, call_next):
        response = await call_next(request)
        # the page and its scripts change with the viewer; browsers check for a new copy on every load
        # instead of reusing a cached one (unchanged files still come back as a quick 304)
        if request.url.path == '/' or request.url.path.startswith('/static/'):
            response.headers['Cache-Control'] = 'no-cache'
        response.headers['Content-Security-Policy'] = CONTENT_SECURITY_POLICY
        response.headers['X-Content-Type-Options'] = 'nosniff'
        response.headers['X-Frame-Options'] = 'DENY'
        response.headers['Referrer-Policy'] = 'no-referrer'
        return response

    app.add_middleware(TrustedHostMiddleware, allowed_hosts=hosts or default_hosts())

    @app.get('/api/config')
    def config():
        config = run.config
        return {
            'run_dir': run.run_dir,
            'wells': run.wells,
            'layouts': {well: run.files.layout(well).mode for well in run.wells},
            'default_well': app.state.default_well,
            'sequencing_channels': list(config['sequencing_channels']),
            'phenotyping_channels': all_phenotyping_channels(config),
            'segmentation_channels': segmentation_channels(config),
            'scale': config['phenotype_scale'] // config['bases_scale'],
        }

    @app.get('/api/well')
    def well(well: str):
        viewer = run.viewer(well)
        return {
            'well': well,
            'tiled': viewer.tables.layout.tiled,
            'seq_cycles': viewer.images.seq_cycles,
            'pt_cycles': viewer.images.pt_cycles,
            'pt_channels': {cycle: viewer.cycle_channels(cycle) for cycle in viewer.images.pt_cycles},
            'tiles': viewer.tables.tiles(),
            'contrast': viewer.contrast(),
        }

    @app.get('/api/dot')
    def dot(well: str, tile: str, index: int, size: int = Query(None, ge=8, le=1024)):
        try:
            return run.viewer(well).dot(tile, index, size)
        except (ValueError, FileNotFoundError) as error:
            raise HTTPException(404, str(error))

    @app.get('/api/cell')
    def cell(well: str, id: int, size: int = Query(None, ge=8, le=1024)):
        try:
            return run.viewer(well).cell(id, size)
        except (ValueError, FileNotFoundError) as error:
            raise HTTPException(404, str(error))

    @app.get('/api/region')
    def region(well: str, tile: str, box: str, cell: int = None):
        try:
            return run.viewer(well).region(tile, box, exclude_cell=cell)
        except (ValueError, FileNotFoundError) as error:
            raise HTTPException(404, str(error))

    @app.get('/api/crop.png')
    def crop(well: str, tile: str, cycle: str, channel: str, box: str, vmin: float = None, vmax: float = None,
            contrast: str = 'global', color: str = Query(None, pattern='^[0-9a-fA-F]{6}$')):
        data = run.viewer(well).crop_png(tile, cycle, channel, box, vmin, vmax, contrast, color)
        return Response(data, media_type='image/png', headers={'Cache-Control': 'private, max-age=3600'})

    @app.get('/api/composite.png')
    def composite(well: str, tile: str, box: str, layer: List[str] = Query([]), contrast: str = 'global'):
        data = run.viewer(well).composite_png(tile, box, layer, contrast)
        return Response(data, media_type='image/png', headers={'Cache-Control': 'private, max-age=3600'})

    @app.get('/')
    def index():
        return FileResponse(os.path.join(STATIC_DIR, 'index.html'))

    app.mount('/static', StaticFiles(directory=STATIC_DIR), name='static')
    return app


def load_token(new=False, path=TOKEN_FILE):
    """ The saved access token, or a new random one saved readable by you only. Keeping it across
    restarts keeps bookmarks and the browser cookie working.
    """
    if not new and os.path.exists(path):
        with open(path) as ifile:
            token = ifile.read().strip()
        if token:
            return token
    token = secrets.token_urlsafe(16)
    os.makedirs(os.path.dirname(path), mode=0o700, exist_ok=True)
    if os.path.exists(path):
        os.remove(path)
    with os.fdopen(os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600), 'w') as ofile:
        ofile.write(token + '\n')
    return token


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('run_dir', help='pipeline run directory (containing config.yaml, rawinput/, stitching/, ...)')
    parser.add_argument('--well', default=None,
            help='well selected when the page opens (default: the first well with all files present)')
    parser.add_argument('--host', default='127.0.0.1',
            help='interface to listen on; keep 127.0.0.1 and tunnel with ssh -J, or 0.0.0.0 if you cannot '
                 'reach compute nodes by ssh (then the token and images cross the cluster network unencrypted)')
    parser.add_argument('--port', type=int, default=8000)
    parser.add_argument('--token', nargs='?', const='new', default=None,
            help='access token; default: the one saved in {} (made on first use). '
                 '"--token new" makes and saves a new one, so old links and cookies stop working'.format(TOKEN_FILE))
    parser.add_argument('--layout', choices=LAYOUTS, default='auto',
            help='tiled (sequencing/<well>_grid<N>/tile*/...), untiled (sequencing/<well>/...), or auto: '
                 'tiled if its tables exist for the well, otherwise untiled if those exist')
    parser.add_argument('--cache-dir', default=None, help='default: <run_dir>/viewer_cache/')
    parser.add_argument('--dots-table', default='cells_quality.csv')
    parser.add_argument('--reads-table', default='cells_reads.csv')
    parser.add_argument('--cells-table', default='cells.csv')
    parser.add_argument('--nuclei-table', default='nuclei.csv', help='nuclei matched to the cell table (optional)')
    parser.add_argument('--prebuild-cache', action='store_true',
            help='convert every tile table to parquet before serving (otherwise done on first view of a tile)')
    args = parser.parse_args()

    # always required: 127.0.0.1 is shared by every user on the node, not just you
    token = args.token if args.token not in (None, 'new') else load_token(new=args.token == 'new')

    run = Run(args.run_dir, cache_dir=args.cache_dir, layout=args.layout, dots_table=args.dots_table,
            reads_table=args.reads_table, cells_table=args.cells_table, nuclei_table=args.nuclei_table)
    default_well = args.well or run.default_well()
    if default_well not in run.wells:
        sys.exit('Well {} not found; this run has: {}'.format(default_well, ', '.join(run.wells) or 'no wells'))
    print('Wells:', ', '.join(run.wells), '| opening', default_well,
            '({})'.format(run.files.layout(default_well).mode), flush=True)

    viewer = None
    try:
        viewer = run.viewer(default_well)
    except MissingFiles as error:
        print(error.message + ':', *[f['path'] for f in error.files if not f['exists']], sep='\n    ', flush=True)

    if args.prebuild_cache and viewer is not None:
        for tile in viewer.tables.tiles():
            start = time.time()
            try:
                viewer.tables._dots_parquet(tile)
                viewer.tables._reads_parquet(tile)
                print('cached', tile, '{:.0f}s'.format(time.time() - start), flush=True)
            except FileNotFoundError as error:
                print('skipping', tile, error, flush=True)

    if viewer is not None:
        viewer.warm_up()
    app = make_app(run, token)
    app.state.default_well = default_well

    import socket
    import uvicorn
    node = socket.gethostname()
    suffix = '/?token=' + token
    print('\nViewer for {} running on {}:{}'.format(run.run_dir, node, args.port))
    if args.host in ('127.0.0.1', 'localhost'):
        print('On your computer run:\n    ssh -N -L {0}:localhost:{0} -J <user>@<login node> <user>@{1}'.format(args.port, node))
    else:
        print('On your computer run:\n    ssh -N -L {0}:{1}:{0} <user>@<login node>'.format(args.port, node))
    print('then open  http://localhost:{}{}\n'.format(args.port, suffix), flush=True)
    uvicorn.run(app, host=args.host, port=args.port, log_level='warning')


if __name__ == '__main__':
    main()
