""" Which wells a run has, and which files the viewer needs to show them.
Wells follow workflow/rules/config.smk: config['wells'] if set (an int n means well1..welln),
otherwise every well detected from the file names under rawinput/ as in config.smk
"""

import os
import re
import glob

import numpy as np


class MissingFiles(Exception):
    """ Raised when files needed to show something don't exist. files is a list of
    {'label', 'path', 'exists'} for everything that was needed, missing or not.
    """
    def __init__(self, message, files):
        super().__init__(message)
        self.message = message
        self.files = files

    def as_dict(self):
        return {
            'detail': self.message,
            'expected': self.files,
            'missing': [f for f in self.files if not f['exists']],
        }


def natural_key(name):
    return [int(part) if part.isdigit() else part for part in re.split(r'(\d+)', name)]


def well_names(config, inputfiles):
    wells = config.get('wells')
    if type(wells) == int:
        wells = ['well{}'.format(i) for i in range(1, wells + 1)]
    if wells is None:
        wells = list(inputfiles)
    return sorted(wells, key=natural_key)


def expected_cycles(inputfiles):
    """ Every cycle found for any well: sequencing cycles in order, then phenotype cycles """
    cycles = set()
    for paths in inputfiles.values():
        cycles.update(paths)
    seq = sorted(c for c in cycles if not c.startswith('P'))
    pt = sorted((c for c in cycles if c.startswith('P')), key=lambda c: (c != 'PT', natural_key(c)))
    return seq + pt


WHOLE_WELL = 'whole'
LAYOUTS = ('auto', 'tiled', 'untiled')


class WellLayout:
    """ Where the tables of a well are, and which part of the well each tile of them covers.
    A tiled run (sequencing/<well>_grid<N>/tileXXxYYy/...) has one tile per grid square. An untiled
    run (sequencing/<well>/..., from targets like sequencing/well1/cells_reads.csv) is treated as a
    single tile, 'whole', covering the whole well. Paths are relative to the run directory.
    """

    def __init__(self, run_dir, config, well, mode='auto', dots_table='cells_quality.csv'):
        if mode not in LAYOUTS:
            raise ValueError('Unknown layout {}; use one of {}'.format(mode, ', '.join(LAYOUTS)))
        self.run_dir = run_dir
        self.config = config
        self.well = well
        self.grid_size = int(config.get('sequencing_grid_size', 5))
        self.grid = '{}_grid{}'.format(well, self.grid_size)
        if mode == 'auto':
            # tiled unless only the untiled tables exist; with neither, the missing files are listed as tiled
            exists = lambda path: os.path.exists(os.path.join(run_dir, path))
            has_tiled = exists(self.grid_positions_file()) and glob.glob(
                    os.path.join(run_dir, self.directory('sequencing'), self.grid, 'tile*', dots_table))
            has_untiled = exists(os.path.join(self.directory('sequencing'), well, dots_table))
            mode = 'untiled' if has_untiled and not has_tiled else 'tiled'
        self.mode = mode
        self.tiled = mode == 'tiled'

    def directory(self, kind):
        return self.config.get(kind + '_dir', kind + '/')

    def tile_dir(self, kind, tile):
        """ Directory of a tile's sequencing or segmentation (kind) tables """
        if self.tiled:
            return os.path.join(self.directory(kind), self.grid, tile)
        return os.path.join(self.directory(kind), self.well)

    def grid_positions_file(self):
        return os.path.join(self.directory('stitching'), self.grid, 'grid_composite.json')

    def positions_file(self):
        """ The file giving the tile sections: the grid composite, or for an untiled run the well
        composite, whose bounds stitch_well crops the whole-well raw.tif to
        """
        if self.tiled:
            return self.grid_positions_file()
        return os.path.join(self.directory('stitching'), self.well, 'composite.json')

    def sections(self):
        """ tile -> (r0, c0, r1, c1) in bases pixels of the well composite """
        from .imaging import load_boxes
        positions, sizes = load_boxes(os.path.join(self.run_dir, self.positions_file()))
        if not self.tiled:
            return {WHOLE_WELL: np.concatenate([positions.min(axis=0), (positions + sizes).max(axis=0)])}
        sections = {}
        for x in range(self.grid_size):
            for y in range(self.grid_size):
                k = x * self.grid_size + y
                sections['tile{:02}x{:02}y'.format(x, y)] = np.concatenate([positions[k], positions[k] + sizes[k]])
        return sections


class RunFiles:
    def __init__(self, run_dir, config, inputfiles, dots_table='cells_quality.csv',
            reads_table='cells_reads.csv', cells_table='cells.csv', nuclei_table='nuclei.csv', layout='auto'):
        # nuclei_table is optional (nuclei outlines are left out without it), so it isn't required here
        self.run_dir = run_dir
        self.config = config
        self.inputfiles = inputfiles
        self.cycles = expected_cycles(inputfiles)
        self.dots_table = dots_table
        self.layout_mode = layout
        self._layouts = {}
        self.tables = {'dots': ('sequencing', dots_table), 'reads': ('sequencing', reads_table), 'cells': ('segmentation', cells_table)}

    def layout(self, well):
        """ The WellLayout of a well, detected once """
        if well not in self._layouts:
            self._layouts[well] = WellLayout(self.run_dir, self.config, well, self.layout_mode, self.dots_table)
        return self._layouts[well]

    def entry(self, label, path):
        return {'label': label, 'path': os.path.relpath(path, self.run_dir) if os.path.isabs(path) else path,
                'exists': os.path.exists(os.path.join(self.run_dir, path))}

    def directory(self, kind):
        return self.config.get(kind + '_dir', kind + '/')

    def well_files(self, well):
        """ Files needed for every image of a well """
        rawinput = self.config.get('rawinput_dir', 'rawinput/')
        stitching = self.directory('stitching')
        paths = self.inputfiles.get(well, {})
        files = [self.entry('pipeline config', 'config.yaml')]
        for cycle in self.cycles:
            if cycle in paths:
                files.append(self.entry('cycle {} images'.format(cycle), paths[cycle]))
            else:
                kind = 'phenotype' if cycle.startswith('P') else 'sequencing'
                files.append({'label': 'cycle {} images'.format(cycle), 'exists': False,
                        'path': '{}… ({} .nd2 for {} not found)'.format(rawinput, kind, well)})
            files.append(self.entry('cycle {} stitching positions'.format(cycle),
                    os.path.join(stitching, well, 'cycle' + cycle, 'composite.json')))
        layout = self.layout(well)
        files.append(self.entry('grid tile positions' if layout.tiled else 'well positions', layout.positions_file()))
        return files

    def tile_file(self, well, tile, kind):
        folder, name = self.tables[kind]
        return os.path.join(self.layout(well).tile_dir(folder, tile), name)

    def tile_files(self, well, tile, kinds):
        labels = {'dots': 'dot table', 'reads': 'cell reads table', 'cells': 'cell table'}
        where = tile if self.layout(well).tiled else well
        return [self.entry('{} {}'.format(where, labels[kind]), self.tile_file(well, tile, kind)) for kind in kinds]

    def require(self, files, message):
        if not all(f['exists'] for f in files):
            raise MissingFiles(message, files)
