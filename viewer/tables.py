""" Table lookups for the viewer.

Coordinate conventions, from the pipeline:
    cells_quality.csv (sequencing/<well>_grid<N>/tileXXxYYy/): position_x is the row and position_y the
        column, in bases pixels relative to that tile's raw.tif. 'cell' is the tile-local cell label,
        which is row cell-1 of that tile's cells.csv (0 = no cell).
    cells.csv (segmentation/<well>_grid<N>/tileXXxYYy/): bbox_x* are rows and bbox_y* columns, in
        phenotype pixels relative to that tile's raw_pt.tif. The index is a well-unique cell id,
        also used as the index of cells_reads.csv.
    cells_reads.csv: dot_indicies is a ':'-joined list of cells_quality.csv indices in that tile.

An untiled run has the same tables in sequencing/<well>/ and segmentation/<well>/, relative to the
whole-well raw.tif / raw_pt.tif, and is handled as one tile ('whole') whose section is the well
composite's bounds (see files.WellLayout). There the cell label is also the cell id.

The large per-tile dot tables are converted once to parquet in the cache directory, split into
row groups so single dots can be fetched without loading the table.
"""

import os
import re
import glob
import base64
import threading
import collections

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq


TILE_RE = re.compile(r'^tile(\d+)x(\d+)y$')


def parse_tile(tile):
    """ Accepts 'tile02x03y', '02x03y' or '2,3' and returns 'tile02x03y', x, y """
    tile = tile.strip()
    match = TILE_RE.match(tile) or TILE_RE.match('tile' + tile)
    if match:
        x, y = int(match.group(1)), int(match.group(2))
    else:
        parts = re.split(r'[,x\s]+', tile.strip('y'))
        if len(parts) != 2 or not all(p.isdigit() for p in parts):
            raise ValueError('Unrecognized tile name ' + repr(tile))
        x, y = map(int, parts)
    return 'tile{:02}x{:02}y'.format(x, y), x, y


MASK_COLUMN_RE = re.compile(r'^mask(\d+(?:\.\d+)?)?$')


def mask_columns(columns):
    """ (scale, column) for every encoded mask column (mask, mask2, mask8, mask16.0, ...),
    finest scale first. A bare 'mask' column is full resolution, scale 1.
    """
    found = []
    for name in columns:
        match = MASK_COLUMN_RE.match(str(name))
        if match:
            found.append((float(match.group(1) or 1), name))
    return sorted(found)


def decode_mask(encoded, size, scale=8):
    """ Decodes a base85 packed binary mask from cells.csv (a mask<scale> column), following
    Cell.decode_mask in starcall/cells.py. size is the bbox size in phenotype pixels.
    """
    data = base64.a85decode(encoded.encode('ascii'))
    arr = np.unpackbits(np.frombuffer(data, np.uint8))
    dims = np.ceil(np.asarray(size) * (1 / scale)).astype(int)
    return arr[:dims.prod()].reshape(dims).astype(bool)


def read_cell_csv(path, columns):
    """ Reads columns of a cells.csv / nuclei.csv plus all of its encoded mask columns """
    header = pd.read_csv(path, nrows=0).columns
    return pd.read_csv(path, index_col=0, usecols=list(columns) + [name for scale, name in mask_columns(header)])


class RunTables:
    def __init__(self, run_dir, well, config, layout, cache_dir=None,
            dots_table='cells_quality.csv', reads_table='cells_reads.csv', cells_table='cells.csv',
            nuclei_table='nuclei.csv'):
        self.run_dir = run_dir
        self.well = well
        self.config = config
        self.scale = config['phenotype_scale'] // config['bases_scale']
        self.layout = layout
        self.grid_size = layout.grid_size
        if layout.tiled and int(config.get('segmentation_grid_size', self.grid_size)) != self.grid_size:
            raise ValueError('segmentation_grid_size and sequencing_grid_size differ; cells_reads tables require them to match')
        self.dots_table, self.reads_table, self.cells_table = dots_table, reads_table, cells_table
        # nuclei.csv is matched to cells.csv: row k is the nucleus of cell k, with the same index
        self.nuclei_table = nuclei_table

        self.cache_dir = cache_dir or os.path.join(run_dir, 'viewer_cache')
        os.makedirs(self.cache_dir, exist_ok=True)

        # tile sections in bases pixels of the well composite, as used by stitch_tile / stitch_well
        self.sections = layout.sections()

        self._locks = collections.defaultdict(threading.Lock)
        self._cells = collections.OrderedDict()
        self._nuclei = collections.OrderedDict()
        self._id_ranges = None

    def tiles(self):
        return sorted(self.sections)

    def sequencing_path(self, tile, name):
        return os.path.join(self.run_dir, self.layout.tile_dir('sequencing', tile), name)

    def segmentation_path(self, tile, name):
        return os.path.join(self.run_dir, self.layout.tile_dir('segmentation', tile), name)

    def section(self, tile, phenotype=False):
        """ (r0, c0, r1, c1) of a grid tile; phenotype sections are scaled like stitch_tile_pt """
        section = self.sections[tile]
        if phenotype:
            # position and size are scaled separately in the pipeline
            pos, size = section[:2] * self.scale, (section[2:] - section[:2]) * self.scale
            section = np.concatenate([pos, pos + size])
        return section

    ##### dots #####

    def _dots_parquet(self, tile):
        """ Converts the tile's dot table to parquet once, keeping only the columns the viewer uses """
        path = os.path.join(self.cache_dir, '{}_{}_{}.parquet'.format(self.well, tile, os.path.splitext(self.dots_table)[0]))
        source = self.sequencing_path(tile, self.dots_table)
        if not os.path.exists(source):
            raise FileNotFoundError(source)
        with self._locks[path]:
            if os.path.exists(path) and os.path.getmtime(path) >= os.path.getmtime(source):
                return path
            header = pd.read_csv(source, nrows=0).columns
            keep = [header[0]] + [c for c in header[1:] if c in ('position_x', 'position_y', 'cell', 'max_seq',
                    'mean_chastity', 'min_chastity') or c.startswith('values_cycle') or c.startswith('chastity_cycle')]
            dtypes = {c: np.float32 for c in keep if c.startswith(('values_', 'chastity_', 'mean_', 'min_', 'position_'))}
            dtypes.update({header[0]: np.int64, 'cell': np.int64, 'max_seq': str})
            tmp = path + '.tmp'
            writer = None
            for chunk in pd.read_csv(source, usecols=keep, dtype=dtypes, chunksize=100000):
                chunk = chunk.rename(columns={header[0]: 'dot_index'})
                table = pa.Table.from_pandas(chunk, preserve_index=False)
                if writer is None:
                    writer = pq.ParquetWriter(tmp, table.schema)
                writer.write_table(table, row_group_size=20000)
            writer.close()
            os.replace(tmp, path)
        return path

    def dots(self, tile, indices):
        """ Rows of the tile's dot table with the given indices, as a DataFrame indexed by dot_index """
        path = self._dots_parquet(tile)
        indices = [int(i) for i in indices]
        table = ds.dataset(path).to_table(filter=ds.field('dot_index').isin(indices))
        return table.to_pandas().set_index('dot_index').reindex(indices).dropna(how='all')

    def dots_in_cell(self, tile, label):
        path = self._dots_parquet(tile)
        table = ds.dataset(path).to_table(filter=ds.field('cell') == int(label))
        return table.to_pandas().set_index('dot_index')

    def dots_in_box(self, tile, box):
        """ Dots of a tile whose position is inside box (r0, c0, r1, c1), bases pixels of the well composite """
        section = self.section(tile)
        r0, c0, r1, c1 = (float(v) for v in (box[0] - section[0], box[1] - section[1], box[2] - section[0], box[3] - section[1]))
        x, y = ds.field('position_x'), ds.field('position_y')
        table = ds.dataset(self._dots_parquet(tile)).to_table(
                columns=['dot_index', 'position_x', 'position_y', 'cell', 'max_seq', 'min_chastity'],
                filter=(x >= r0) & (x < r1) & (y >= c0) & (y < c1))
        return table.to_pandas().set_index('dot_index')

    def num_cycles(self, dots):
        return sorted({c[len('chastity_cycle'):] for c in dots.columns if c.startswith('chastity_cycle')})

    ##### cells #####

    def _cell_table(self, tile):
        with self._locks['cells' + tile]:
            if tile in self._cells:
                self._cells.move_to_end(tile)
                return self._cells[tile]
            table = read_cell_csv(self.segmentation_path(tile, self.cells_table),
                    ['Unnamed: 0', 'bbox_x1', 'bbox_y1', 'bbox_x2', 'bbox_y2'])
            self._cells[tile] = table
            while len(self._cells) > 6:
                self._cells.popitem(last=False)
        return table

    def _cell_id_ranges(self):
        """ (min id, max id, tile) for every tile, reading only the index column of each cells.csv """
        if self._id_ranges is None:
            ranges = []
            for tile in self.tiles():
                path = self.segmentation_path(tile, self.cells_table)
                if not os.path.exists(path): continue
                index = pd.read_csv(path, usecols=[0]).iloc[:, 0]
                if len(index):
                    ranges.append((int(index.min()), int(index.max()), tile))
            self._id_ranges = ranges
        return self._id_ranges

    def cell_by_id(self, cell_id):
        """ Returns (tile, label, row) for a well-unique cell id """
        cell_id = int(cell_id)
        for low, high, tile in self._cell_id_ranges():
            if low <= cell_id <= high:
                table = self._cell_table(tile)
                if cell_id in table.index:
                    label = table.index.get_loc(cell_id) + 1
                    return tile, label, table.loc[cell_id]
        raise KeyError('Cell id {} not found in any {}'.format(cell_id, self.cells_table))

    def cell_by_label(self, tile, label):
        """ Returns (cell id, row) for a tile-local cell label from the dot table """
        table = self._cell_table(tile)
        return int(table.index[int(label) - 1]), table.iloc[int(label) - 1]

    def has_nuclei(self, tile):
        return os.path.exists(self.segmentation_path(tile, self.nuclei_table))

    def nucleus_table(self, tile):
        """ The tile's nuclei.csv (bbox, mask columns and orig_index, the label in the unmatched nuclei mask image) """
        with self._locks['nuclei' + tile]:
            if tile in self._nuclei:
                self._nuclei.move_to_end(tile)
                return self._nuclei[tile]
            table = read_cell_csv(self.segmentation_path(tile, self.nuclei_table),
                    ['Unnamed: 0', 'orig_index', 'bbox_x1', 'bbox_y1', 'bbox_x2', 'bbox_y2'])
            self._nuclei[tile] = table
            while len(self._nuclei) > 6:
                self._nuclei.popitem(last=False)
        return table

    def nucleus_by_label(self, tile, label):
        """ Row of the nucleus matched to a tile-local cell label, or None """
        if not self.has_nuclei(tile):
            return None
        table = self.nucleus_table(tile)
        return table.iloc[int(label) - 1] if 0 < int(label) <= len(table) else None

    def cells_in_box(self, tile, box_pt, kind='cells'):
        """ (labels, rows) of the tile's cells (or their nuclei, kind='nuclei') whose bbox overlaps
        box_pt (phenotype pixels of the well composite). Labels are tile-local cell labels.
        """
        table = self._cell_table(tile) if kind == 'cells' else self.nucleus_table(tile)
        section = self.section(tile, phenotype=True)
        r0, c0, r1, c1 = box_pt[0] - section[0], box_pt[1] - section[1], box_pt[2] - section[0], box_pt[3] - section[1]
        hits = np.nonzero(((table.bbox_x1 < r1) & (table.bbox_x2 > r0) & (table.bbox_y1 < c1) & (table.bbox_y2 > c0)).to_numpy())[0]
        return hits + 1, table.iloc[hits]

    def cell_box(self, tile, row):
        """ Cell bbox (r0, c0, r1, c1) in phenotype pixels of the well composite """
        section = self.section(tile, phenotype=True)
        return np.array([row.bbox_x1, row.bbox_y1, row.bbox_x2, row.bbox_y2]) + np.tile(section[:2], 2)

    def cell_mask(self, row):
        """ Full resolution (phenotype pixels) boolean mask of a cell over its bbox, from the finest
        mask column the row has, or None
        """
        size = (int(row.bbox_x2 - row.bbox_x1), int(row.bbox_y2 - row.bbox_y1))
        for scale, name in mask_columns(row.index):
            encoded = row[name]
            if isinstance(encoded, str) and encoded:
                small = decode_mask(encoded, size, scale)
                # nearest neighbour upscaling, pixel i of the bbox is pixel floor(i / scale) of the mask
                rows = np.minimum((np.arange(size[0]) / scale).astype(int), small.shape[0] - 1)
                cols = np.minimum((np.arange(size[1]) / scale).astype(int), small.shape[1] - 1)
                return small[np.ix_(rows, cols)]
        return None

    ##### reads #####

    def _reads_parquet(self, tile):
        path = os.path.join(self.cache_dir, '{}_{}_{}.parquet'.format(self.well, tile, os.path.splitext(self.reads_table)[0]))
        source = self.sequencing_path(tile, self.reads_table)
        if not os.path.exists(source):
            raise FileNotFoundError(source)
        with self._locks[path]:
            if not (os.path.exists(path) and os.path.getmtime(path) >= os.path.getmtime(source)):
                table = pd.read_csv(source, index_col=0, dtype=str)
                table.index = table.index.astype(np.int64)
                table.index.name = 'cell_id'
                table.to_parquet(path + '.tmp', row_group_size=20000)
                os.replace(path + '.tmp', path)
        return path

    def cell_reads(self, tile, cell_id):
        """ The cells_reads row of a cell as a dict of non-empty string values, or {} """
        table = ds.dataset(self._reads_parquet(tile)).to_table(filter=ds.field('cell_id') == int(cell_id)).to_pandas()
        if len(table) == 0:
            return {}
        row = table.iloc[0].drop(labels=['cell_id'], errors='ignore')
        return {k: v for k, v in row.items() if isinstance(v, str) and v != ''}
