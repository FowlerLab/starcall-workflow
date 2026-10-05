""" Lazy image access for the viewer.

Crops are read straight out of the raw nd2 files using the solved per-cycle composites
in stitching/, so raw.tif / raw_pt.tif are not needed. Each cycle's composite gives an
integer (row, col) origin for every nd2 frame (frame j == box j); tiles are never resized
or flipped, so a crop is just a paste of the few frames that overlap it.

Cycles listed as rotated in stitching/<well>/rotation.csv (detect_rotation in
workflow/rules/alignment.smk) had every frame rotated about its centre by applied_deg before
alignment and stitching, keeping its shape, so the composite positions are those of the rotated
frames. The same rotation is applied here to the part of each frame that is pasted, matching
starcall.rotation.rotate_frame. Without rotation.csv (runs from before rotation correction)
no frame is rotated.

The paste reproduces how stitch_well_section (workflow/rules/stitching.smk) builds raw.tif:
tiles are clipped to the grid tile section, then merged with constitch's EfficientNearestMerger
(constitch/merging.py), where the pixel farthest from its (clipped) tile edge wins, distances are
capped at 255, and ties go to tiles added in later stitching batches (constitch/composite.py stitch()).
"""

import os
import json
import threading
import collections

import numpy as np
import yaml


def load_config(run_dir):
    """ Loads config.yaml of a pipeline run, filling missing keys from default-config.yaml """
    config = {}
    for name in ('default-config.yaml', 'config.yaml'):
        path = os.path.join(run_dir, name)
        if os.path.exists(path):
            with open(path) as ifile:
                config.update(yaml.safe_load(ifile) or {})
    return config


def find_input_files(run_dir, config):
    """ Maps well -> cycle -> nd2/tif path, with the same rules as workflow/rules/config.smk
    (sorted file paths under rawinput/, phenotype_date in the path marks phenotype cycles).
    """
    rawinput_dir = os.path.join(run_dir, config.get('rawinput_dir', 'rawinput/'))
    phenotype_date = config.get('phenotype_date', 'phenotype')

    cycles = config.get('cycles')
    if type(cycles) == int:
        cycles = ['{:02}'.format(i) for i in range(cycles)]
    elif cycles is not None and type(cycles[0]) == int:
        cycles = ['{:02}'.format(i) for i in cycles]

    pt_cycles = config.get('phenotype_cycles')
    wells = config.get('wells')
    if type(wells) == int:
        wells = ['well{}'.format(i) for i in range(1, wells + 1)]
    if type(pt_cycles) == int:
        pt_cycles = (['PT'] + ['P{}'.format(i) for i in range(1, pt_cycles)])[:pt_cycles]

    paths = []
    for root, dirs, files in os.walk(rawinput_dir, followlinks=True):
        paths.extend(root + '/' + filename for filename in files)
    paths = sorted(paths)

    inputfiles = {}
    for path in paths:
        if not (path.endswith('.tif') or path.endswith('.tiff') or path.endswith('.nd2')):
            continue
        # paths relative to rawinput/, so run directory names can't match well or phenotype patterns
        relpath = path[len(rawinput_dir):]
        if wells is not None:
            matching = [well for well in wells if any(option in relpath for option in
                    [well, well[:4].replace('well', 'Well') + well[4:], 'well' + well, 'Well' + well])]
            if len(matching) == 0: continue
            well = sorted(matching, key=len)[-1]
        else:
            if not (relpath.count('well') or relpath.count('Well')): continue
            well = relpath.replace('Well', 'well')
            well = well[well.index('well'):].split('_')[0]
        inputfiles.setdefault(well, []).append((relpath, path))

    result = {}
    for well, wellpaths in inputfiles.items():
        cyclepaths = {}
        index = pt_index = 0
        for relpath, path in wellpaths:
            if phenotype_date in relpath:
                if pt_cycles is not None:
                    if pt_index >= len(pt_cycles): continue
                    cycle = pt_cycles[pt_index]
                else:
                    cycle = 'PT' if pt_index == 0 else 'P{}'.format(pt_index)
                pt_index += 1
            else:
                if cycles is not None:
                    if index >= len(cycles): continue
                    cycle = cycles[index]
                else:
                    cycle = '{:02}'.format(index)
                index += 1
            cyclepaths[cycle] = path
        result[well] = cyclepaths
    return result


def load_boxes(path):
    """ Reads a constitch composite json into (positions, sizes) int arrays of shape (N, 2) """
    with open(path) as ifile:
        positions, sizes = json.load(ifile)['composite']['boxes']
    positions = np.array(positions, dtype=np.int64)[:, :2]
    sizes = np.array(sizes, dtype=np.int64)[:, :2]
    return positions, sizes


def load_rotations(path):
    """ cycle -> rotation (degrees) applied to the frames of the cycle, from the rotation.csv
    written by detect_rotation. Empty if the run has no rotation.csv. """
    if not os.path.exists(path):
        return {}
    import csv
    with open(path, newline='') as ifile:
        return {row['cycle']: float(row['applied_deg']) for row in csv.DictReader(ifile)}


def rotate_region(frame, angle, rows, cols, order=1):
    """ Rows rows[0]:rows[1] and columns cols[0]:cols[1] of rotate_frame(frame, angle, order)
    (starcall/rotation.py, scipy.ndimage.rotate about the frame centre with mode='nearest'),
    for a (C, H, W) frame, without rotating the rest of the frame """
    if angle == 0:
        return frame[:, rows[0]:rows[1], cols[0]:cols[1]]
    import scipy.ndimage
    import scipy.special
    # the transform scipy.ndimage.rotate uses, with the output shifted to the region's origin
    c, s = scipy.special.cosdg(angle), scipy.special.sindg(angle)
    matrix = np.array([[c, s], [-s, c]])
    center = (np.array(frame.shape[1:]) - 1) / 2
    offset = center - matrix @ center + matrix @ np.array([rows[0], cols[0]])
    shape = (rows[1] - rows[0], cols[1] - cols[0])
    out = np.empty((frame.shape[0],) + shape, frame.dtype)
    for k in range(frame.shape[0]):
        scipy.ndimage.affine_transform(frame[k], matrix, offset, shape, out[k], order=order, mode='nearest')
    return out


def stitch_batches(positions, sizes):
    """ Order in which constitch's stitch() adds tiles: it repeatedly takes, in index order,
    every tile not overlapping one already taken in the current pass. Tiles in later passes
    win ties in the nearest merger, so the pass number is returned for every tile.
    """
    point1, point2 = positions, positions + sizes
    batch = np.full(len(positions), -1)
    left = list(range(len(positions)))
    cur = 0
    while left:
        skipped, taken = [], []
        for i in left:
            if taken:
                t = np.array(taken)
                overlaps = np.all((point1[t] < point2[i]) & (point1[i] < point2[t]), axis=1)
                if overlaps.any():
                    skipped.append(i)
                    continue
            taken.append(i)
            batch[i] = cur
        left = skipped
        cur += 1
    return batch


class CycleImages:
    """ One cycle of one well: the nd2 file plus its composite positions, and the rotation
    (degrees) applied to its frames, see detect_rotation """

    def __init__(self, path, composite_path, angle=0, rotation_order=1):
        self.path = path
        self.angle = angle
        self.rotation_order = rotation_order
        self.positions, self.sizes = load_boxes(composite_path)
        self.batches = stitch_batches(self.positions, self.sizes)
        self._file = None
        self.lock = threading.Lock()

    def file(self):
        if self._file is None:
            import nd2
            self._file = nd2.ND2File(self.path)
        return self._file

    def read_frame(self, index):
        """ Frame index as (channels, rows, cols); nd2 Y is the row axis, as in the pipeline """
        if self.path.endswith('.nd2'):
            return self.file().read_frame(int(index))
        import tifffile
        if self._file is None:
            self._file = tifffile.memmap(self.path, mode='r')
        return self._file[index]

    def read_region(self, index, rows, cols):
        """ Rows rows[0]:rows[1] and columns cols[0]:cols[1] of frame index, as stitched:
        rotated by the cycle's rotation """
        return rotate_region(self.read_frame(index), self.angle, rows, cols, self.rotation_order)

    def crop(self, box, section=None):
        """ Crops box = (r0, c0, r1, c1) (in this cycle's composite pixels) from all channels.
        section = (r0, c0, r1, c1) is the grid tile the matching raw.tif was stitched from;
        tiles are clipped to it before measuring edge distances, exactly as in the pipeline.
        Pixels of the box outside the section use unclipped tiles. Pixels with no tile are NaN.
        """
        r0, c0, r1, c1 = box
        shape = (r1 - r0, c1 - c0)
        best = np.full(shape, -1, dtype=np.int64)
        image = None

        regions = [((r0, c0, r1, c1), None)]
        if section is not None:
            inner = (max(r0, section[0]), max(c0, section[1]), min(r1, section[2]), min(c1, section[3]))
            if inner[0] < inner[2] and inner[1] < inner[3]:
                regions.append((inner, section))

        with self.lock:
            for (a0, b0, a1, b1), clip in regions:
                point1, point2 = self.positions, self.positions + self.sizes
                if clip is not None:
                    point1 = np.maximum(point1, clip[:2])
                    point2 = np.minimum(point2, clip[2:])
                hits = np.nonzero(np.all((point1 < (a1, b1)) & ((a0, b0) < point2), axis=1))[0]
                if clip is not None:
                    best[a0-r0:a1-r0, b0-c0:b1-c0] = -1
                    if image is not None:
                        image[:, a0-r0:a1-r0, b0-c0:b1-c0] = np.nan

                for j in hits:
                    p1, p2 = point1[j], point2[j]
                    q0, q1 = max(a0, p1[0]), min(a1, p2[0])
                    s0, s1 = max(b0, p1[1]), min(b1, p2[1])
                    rows, cols = np.arange(q0, q1), np.arange(s0, s1)
                    rdist = np.minimum(np.minimum(rows - p1[0], p2[0] - 1 - rows), 254)
                    cdist = np.minimum(np.minimum(cols - p1[1], p2[1] - 1 - cols), 254)
                    dist = np.minimum(rdist[:, None], cdist[None, :]) + 1
                    # ties go to later batches; tiles in one batch never overlap
                    key = dist * 1000 + self.batches[j]

                    origin = self.positions[j]
                    tile = self.read_region(j, (q0 - origin[0], q1 - origin[0]), (s0 - origin[1], s1 - origin[1]))
                    if image is None:
                        image = np.full((tile.shape[0],) + shape, np.nan, dtype=np.float32)
                    dest = (slice(q0 - r0, q1 - r0), slice(s0 - c0, s1 - c0))
                    mask = key > best[dest]
                    image[(slice(None),) + dest][:, mask] = tile[:, mask]
                    best[dest][mask] = key[mask]

        if image is None:
            num_channels = self.read_frame(0).shape[0]
            image = np.full((num_channels,) + shape, np.nan, dtype=np.float32)
        return image

    def coverage(self, box):
        """ Fraction of box (r0, c0, r1, c1) that some frame of this cycle covers. Parts of the well
        outside every frame were not imaged in this cycle (NaN in raw.tif).
        """
        r0, c0, r1, c1 = box
        covered = np.zeros((r1 - r0, c1 - c0), bool)
        point1, point2 = self.positions, self.positions + self.sizes
        for j in np.nonzero(np.all((point1 < (r1, c1)) & ((r0, c0) < point2), axis=1))[0]:
            covered[max(point1[j, 0], r0) - r0:min(point2[j, 0], r1) - r0,
                    max(point1[j, 1], c0) - c0:min(point2[j, 1], c1) - c0] = True
        return float(covered.mean())

    def sample_pixels(self, num_frames=3, step=8):
        """ Subsampled pixels from a few frames spread across the well, for default contrast """
        indices = np.linspace(0, len(self.positions) - 1, num_frames + 2).astype(int)[1:-1]
        with self.lock:
            frames = [self.read_frame(j)[:, ::step, ::step] for j in indices]
        return np.concatenate([frame.reshape(frame.shape[0], -1) for frame in frames], axis=1)

    def close(self):
        if self._file is not None and hasattr(self._file, 'close'):
            self._file.close()
        self._file = None


class RunImages:
    """ All cycles of one well of a pipeline run """

    def __init__(self, run_dir, well, config, paths=None, crop_cache_bytes=512 * 2**20):
        self.run_dir = run_dir
        self.well = well
        self.config = config
        self.scale = config['phenotype_scale'] // config['bases_scale']
        stitching_dir = os.path.join(run_dir, config.get('stitching_dir', 'stitching/'))

        if paths is None:
            paths = find_input_files(run_dir, config).get(well)
        if not paths:
            raise ValueError('No input files found for {} under {}'.format(well, config.get('rawinput_dir', 'rawinput/')))
        self.rotations = load_rotations(os.path.join(stitching_dir, well, 'rotation.csv'))
        rotation_order = (config.get('stitching') or {}).get('rotation', {}).get('interpolation_order', 1)
        self.cycles = {}
        for cycle, path in paths.items():
            composite_path = os.path.join(stitching_dir, well, 'cycle' + cycle, 'composite.json')
            if not os.path.exists(composite_path):
                raise FileNotFoundError('Missing solved composite ' + composite_path)
            self.cycles[cycle] = CycleImages(path, composite_path, self.rotations.get(cycle, 0), rotation_order)

        self.seq_cycles = sorted(c for c in self.cycles if not c.startswith('P'))
        self.pt_cycles = [c for c in self.cycles if c.startswith('P')]

        self._cache = collections.OrderedDict()
        # limited by memory rather than count: one large phenotype crop is hundreds of MB
        self._cache_bytes = crop_cache_bytes
        self._cache_used = 0
        self._cache_lock = threading.Lock()

    def is_pt(self, cycle):
        return cycle in self.pt_cycles

    def crop(self, cycle, box, section=None):
        """ All channels of box (r0, c0, r1, c1) of a cycle, in that cycle's pixels, cached """
        key = (cycle, tuple(box), None if section is None else tuple(section))
        with self._cache_lock:
            if key in self._cache:
                self._cache.move_to_end(key)
                return self._cache[key]
        image = self.cycles[cycle].crop(box, section)
        with self._cache_lock:
            if key in self._cache:
                self._cache_used -= self._cache.pop(key).nbytes
            self._cache[key] = image
            self._cache_used += image.nbytes
            while self._cache and self._cache_used > self._cache_bytes:
                self._cache_used -= self._cache.popitem(last=False)[1].nbytes
        return image

    def close(self):
        for cycle in self.cycles.values():
            cycle.close()
