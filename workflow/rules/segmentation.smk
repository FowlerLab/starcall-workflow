import os
import glob
import re

wildcard_constraints:
    output_dir = '|'.join([sequencing_dir, segmentation_dir, phenotyping_dir]),
    path_nogrid2 = '((?!_grid)[^.])*',
    grid = '|_grid{}'.format(config.get('segmentation_grid_size', 1)),
    unmatched = '_unmatched' if config['segmentation'].get('match_masks', False) else '',

rule segment_nuclei:
    """ Uses Stardist to segment the nuclei of cells in the phenotyping
    images. Takes one phenotyping channel as a nuclear channel. Outputs
    an image with integer masks for each cell.

    Params:
        nuclearchannel: The phenotyping channel to perform segmentation on, should be
            an integer index or one of the phenotyping channels specified in config.yaml
    """
    input:
        (segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/corrected_pt.tif'
                if config['segmentation']['use_corrected'] else
                segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/raw_pt.tif'),
    output:
        # grid is included twice as segmentation on grid tiles needs to be merged, so the output is marked with '_grid'
        segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/nuclei{nuclearchannel}_mask{unmatched}{grid}.tif',
    params:
        nuclearchannel = parse_param('nuclearchannel', config['segmentation']['channels'][0]),
        method = config['segmentation']['nuclei_method'],
    wildcard_constraints:
        nuclearchannel = '|_nuclearchannel' + phenotyping_channel_regex,
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 64 + 10000,
        #cuda = 1,
    threads: 2
    run:
        import numpy as np
        import tifffile
        import starcall.segmentation
        import skimage.segmentation

        nuclearchannel = channel_index_phenotyping(params.nuclearchannel)

        data = tifffile.memmap(input[0], mode='r')
        if data.shape[3] < 32:
            data = data.transpose(3,0,1,2)
        #data = data.reshape(-1, *data.shape[2:])

        dapi = data[nuclearchannel[0],nuclearchannel[1]]
        if np.all(dapi == 0):
            tifffile.imwrite(output[0], data[0])
        else:
            del data
            nuclei = starcall.segmentation.segment_nuclei(dapi, method=params.method)
            #nuclei, fmap, rmap = skimage.segmentation.relabel_sequential(skimage.segmentation.clear_border(nuclei))
            nuclei = starcall.segmentation.filter_segmentation(nuclei)
            debug ('Found', nuclei.max(), 'nuclei')
            tifffile.imwrite(output[0], nuclei)


rule segment_cells:
    """ Segments phenotype images to find cell boundaries, used to group reads and calculate cell phenotype features.
    Takes two channels as input, a nuclear channel and a cytoplasm channel.

    Params:
        nuclearchannel, cytochannel: The phenotyping channels to perform segmentation on, should be
            an integer index or one of the phenotyping channels specified in config.yaml
    """
    input:
        (segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/corrected_pt.tif'
                if config['segmentation']['use_corrected'] else
                segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/raw_pt.tif'),
    output:
        segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/cells{diameter}{nuclearchannel}{cytochannel}_mask{unmatched}{grid}.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 20 + 10000,
        #cuda = 1,
    params:
        diameter = parse_param('diameter', config['segmentation']['diameter']),
        nuclearchannel = parse_param('nuclearchannel', config['segmentation']['channels'][0]),
        cytochannel = parse_param('cytochannel', config['segmentation']['channels'][1]),
        method = config['segmentation']['cells_method'],
    wildcard_constraints:
        diameter = '|_diameter\d+',
        nuclearchannel = '|_nuclearchannel' + phenotyping_channel_regex,
        cytochannel = '|_cytochannel' + phenotyping_channel_regex,
    threads: 2
    run:
        import numpy as np
        import starcall.segmentation
        import tifffile
        import logging
        import skimage.segmentation

        nuclearchannel = channel_index_phenotyping(params.nuclearchannel)
        cytochannel = channel_index_phenotyping(params.cytochannel)

        use_gpu = hasattr(resources, 'cuda') and resources.cuda == 1

        debug ('channel indices', nuclearchannel, cytochannel)

        logging.basicConfig(level=logging.INFO)

        data = tifffile.memmap(input[0], mode='r')
        debug (data.shape)
        if data.shape[3] < 32:
            data = data.transpose(3,0,1,2)
        debug (data.shape)
        #data = data.reshape(-1, *data.shape[2:])

        debug(data.shape)

        dapi = data[nuclearchannel[0],nuclearchannel[1]]
        cyto = data[cytochannel[0],cytochannel[1]]
        if np.all(dapi == 0) or np.all(cyto == 0):
            tifffile.imwrite(output[0], data[0])
        else:
            del data
            #del full_well

            cells = starcall.segmentation.segment_cells(
                cyto, dapi,
                method = params.method,
                diameter = params.diameter,
                gpu = use_gpu,
            )
            #cells, fmap, rmap = skimage.segmentation.relabel_sequential(skimage.segmentation.clear_border(cells))
            cells = starcall.segmentation.filter_segmentation(cells)

            debug ('Found', cells.max(), 'cells')

            tifffile.imwrite(output[0], cells)#, compression='deflate')


rule expand_segmentation:
    """ Expands the segmentation by a specified pixel amount, performed with skimage.segmentation.expand_labels
    """
    input:
        #cells = segmentation_dir + '{path}/{segmentation_type}_mask_unmerged.tif',
        cells = segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/{segmentation_type}_mask{downscaled}{unmatched}{grid}.tif',
    output:
        #cells = segmentation_dir + '{path}/{segmentation_type}expanded{size,\d+}_mask_unmerged.tif',
        cells = segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/{segmentation_type}expanded{size,\d+}_mask{downscaled}{unmatched}{grid}.tif',
    wildcard_constraints:
        downscaled = '|_downscaled',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 2 + 5000,
    run:
        import tifffile
        import numpy as np
        import skimage.segmentation

        size = int(wildcards.size)

        cells = tifffile.imread(input.cells)
        cells = skimage.segmentation.expand_labels(cells, size)
        tifffile.imwrite(output.cells, cells)


rule segment_cells_bases:
    """ Segments sequencing images to find cell boundaries, used to group reads and calculate cell phenotype features.
    Takes one channel as input, a nuclear channel. The cytoplasm channel is calculated using starcall.segmentation.estimate_cyto

    Params:
        nuclearchannel: The sequencing channel to perform segmentation on, should be
            an integer index or one of the sequencing channels specified in config.yaml
    """
    input:
        segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/raw.tif'
        #segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/cycle' + cycles[-1] + '/raw.tif'
    output:
        segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/cellsbases{diameter}{nuclearchannel}_mask_downscaled{unmatched}{grid}.tif',
    params:
        diameter = parse_param('diameter', config['segmentation']['diameter']),
        nuclearchannel = parse_param('nuclearchannel', config['segmentation']['channels'][0]),
        method = config['segmentation']['cells_method'],
    wildcard_constraints:
        diameter = '|_diameter\d+',
        nuclearchannel = '|_nuclearchannel' + sequencing_channel_regex,
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 16 + 10000,
        #cuda = 1,
    threads: 8
    run:
        import numpy as np
        import starcall.segmentation
        import tifffile
        import skimage.segmentation

        diameter = params.diameter
        nuclearchannel = channel_index(params.nuclearchannel, kind='sequencing')

        use_gpu = hasattr(resources, 'cuda') and resources.cuda == 1

        full_well = tifffile.memmap(input[0], mode='r')
        data = full_well[cycles.index(cellpose_cycle)].astype(np.float32)

        if np.all(data == 0):
            tifffile.imwrite(output[0], data[0])
            tifffile.imwrite(output[1], data[0])
        else:
            dapi = data[nuclearchannel]
            #cyto = data[2]
            cyto = starcall.segmentation.estimate_cyto(data[sequencing_channels_slice])
            del data
            del full_well

            cells = starcall.segmentation.segment_cells(
                cyto, dapi,
                method = params.method,
                diameter = config['segmentation']['diameter'] * bases_scale // phenotype_scale,
                gpu=use_gpu,
            )
            #cells, fmap, rmap = skimage.segmentation.relabel_sequential(skimage.segmentation.clear_border(cells))
            cells = starcall.segmentation.filter_segmentation(cells)

            debug(f'found {cells.max()} cells ')
            tifffile.imwrite(output[0], cells)#, compression='deflate')


rule segment_nuclei_bases:
    """ Uses Stardist to segment the nuclei of cells from the sequencing
    images. Takes one channel as a nuclear channel. Outputs
    an image with integer masks for each cell.

    Params:
        nuclearchannel: The sequencing channel to perform segmentation on, should be
            an integer index or one of the sequencing channels specified in config.yaml
    """
    input:
        segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/raw.tif'
    output:
        segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/nucleibases{nuclearchannel}_mask_downscaled{unmatched}{grid}.tif',
    params:
        nuclearchannel = parse_param('nuclearchannel', config['segmentation']['channels'][0])
    wildcard_constraints:
        nuclearchannel = '|_nuclearchannel' + sequencing_channel_regex,
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 16 + 10000
        #cuda=1
    threads: 8
    run:
        import numpy as np
        import starcall.segmentation
        import tifffile
        import skimage.segmentation

        nuclearchannel = channel_index(params.nuclearchannel, kind='sequencing')

        full_well = tifffile.memmap(input[0], mode='r')
        data = full_well[-1]

        if np.all(data == 0):
            tifffile.imwrite(output[0], data[0])
            tifffile.imwrite(output[1], data[0])
        else:
            dapi = data[nuclearchannel]
            del data
            del full_well

            nuclei = starcall.segmentation.segment_nuclei(dapi)
            #nuclei, fmap, rmap = skimage.segmentation.relabel_sequential(skimage.segmentation.clear_border(nuclei))
            nuclei = starcall.segmentation.filter_segmentation(nuclei)
            debug(f'found {nuclei.max()} nuclei ')

            tifffile.imwrite(output[0], nuclei)#, compression='deflate')


rule downscale_segmentation:
    """ Rescales cell segmentation to match the scale of the sequencing images. The
    scale ratio is specified in config.yaml, with phenotype_scale and sequencing_scale.
    """
    input:
        segmentation_dir + '{path}/{segmentation_type}_mask.tif',
    output:
        segmentation_dir + '{path}/{segmentation_type}_mask_downscaled.tif',
    run:
        import tifffile
        import skimage.transform

        mask = tifffile.imread(input[0])

        if wildcards.segmentation_type.count('bases') != 0:
            tifffile.imwrite(output[0], mask)
        else:
            tifffile.imwrite(output[0], skimage.transform.rescale(mask, bases_scale/phenotype_scale, order=0))



rule tabulate_cells:
    """ Simple information is recorded about the segmented cells, such as position, bbox.
    Also stores masks of each cell, downscaled so they take up a reasonable amount of space.
    These binary masks are packed 8 bits per byte then encoded using base85 to allow them to be stored
    in a csv file. The encoding and decoding of these masks is handled with the pandas extention 'cells'
    defined in starcall.cells
    If the input segmentation is sequencing image scale (contains '_downscaled' in the path), the bboxes
    and positions of cells is converted to phenotype scale, to ensure all cell tables are consistently
    phenotype scale.
    """
    input:
        cells = lambda wildcards: (
                segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/{segmentation_type}_mask{unmatched}{grid}.tif'
                if 'bases' not in wildcards.segmentation_type else
                segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/{segmentation_type}_mask_downscaled{unmatched}{grid}.tif'),
        #cells = segmentation_dir + '{path}/{segmentation_type}_mask.tif',
    output:
        #table = segmentation_dir + '{path}/{segmentation_type}{unmerged}.csv',
        table = segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/{segmentation_type}{unmatched}{grid}.csv',
        #table = segmentation_dir + '{path}/{segmentation_type}.csv',
    #wildcard_constraints:
        #unmerged = '_unmatched(|_grid\d+)' if config['segmentation'].get('match_masks', False) else '(|_grid\d+)',
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 1.5
    run:
        import numpy as np
        import tifffile
        import skimage.measure
        import starcall.sequencing
        import starcall.cells
        import pandas

        cells = tifffile.imread(input.cells)
        table = starcall.cells.make_cell_table(cells)
        masks_scale = 8

        if 'bases' in wildcards.segmentation_type:
            newscale = (masks_scale * phenotype_scale / bases_scale)
            #masks_scale = round(newscale * bases_scale / phenotype_scale)
            debug (masks_scale, newscale)

            table.cells.rescale_masks(masks_scale)
            table.cells.bboxes[:,:2] = np.ceil(table.cells.bboxes[:,:2] * phenotype_scale / bases_scale).astype(int)
            table.cells.bboxes[:,2:] = np.floor(table.cells.bboxes[:,2:] * phenotype_scale / bases_scale).astype(int)
            #table.cells.bboxes[:] *= phenotype_scale
            #table.cells.bboxes[:] //= bases_scale

            table['mask{}'.format(newscale)] = table['mask{}'.format(masks_scale)]
            table = table.drop('mask{}'.format(masks_scale), axis=1)

        else:
            table.cells.rescale_masks(masks_scale)

        table.to_csv(output.table)


rule plot_cells:
    """ Plots cells using a cell table. If the cell table has cell masks saved by tabulate_cells
    the masks are used, otherwise the bounding boxes of each cell is used.
    """
    input:
        table = segmentation_dir + '{path}/{segmentation_type}{unmerged}.csv',
    output:
        plot = qc_dir + '{path}/{segmentation_type}{unmerged,(|_unmatched)(|_grid\d*)}.svg',
    run:
        import pandas
        import numpy as np
        import matplotlib.pyplot as plt
        import starcall.cells

        table = pandas.read_csv(input.table, index_col=0)

        fig, axes = plt.subplots(figsize=(25, 25))

        table.cells.plot(axes, masks=True)

        #for i in table.index:
            #x1, y1, x2, y2 = table['bbox_x1'][i], table['bbox_y1'][i], table['bbox_x2'][i], table['bbox_y2'][i]
            #axes.plot([x1, x1, x2, x2, x1], [y1, y2, y2, y1, y1], color='C{}'.format(i%10))

        fig.savefig(output.plot)


def neighboring_tables(wildcards):
    x, y = int(wildcards.x), int(wildcards.y)
    grid_size = int(wildcards.grid_size)
    tiles = [(i, j) for i in range(x-1, x+1) for j in range(y-1, y+2)
            if (i * grid_size + j) < (x * grid_size + y) and i >= 0 and 0 <= j < grid_size]
    return [segmentation_dir + '{well}_grid{grid_size}/'
                 + 'tile{:02}x{:02}y'.format(i, j)
                 + '/{segmentation_type}{unmatched}.csv' for i, j in tiles]

rule drop_duplicate_cells:
    """ Removes duplicate cells in the overlapping region between neighboring cells. Cells that are closer to the center of
    another tile in the grid are removed. Afterwards cells that overlap with cells from other tiles are removed as well. Overlap
    between cells is only checked for tiles with an index less than the current tile.
    """
    input:
        table = segmentation_dir + '{well}_grid{grid_size}/tile{x}x{y}y/{segmentation_type}{unmatched}_grid{grid_size}.csv',
        edge_tables = neighboring_tables,
        composite = stitching_dir + '{well}_grid{grid_size}/grid_composite.json',
    output:
        table = segmentation_dir + '{well}_grid{grid_size}/tile{x}x{y}y/{segmentation_type}{unmatched}.csv',
    run:
        import numpy as np
        import constitch
        import pandas
        import sklearn.neighbors
        import starcall.cells

        overlap_threshold = 0.5
        x, y, grid_size = int(wildcards.x), int(wildcards.y), int(wildcards.grid_size)
        index = x * grid_size + y

        composite = constitch.load(input.composite)

        # scaling from base images to phenotype
        composite.boxes.positions[:,:2] *= phenotype_scale
        composite.boxes.positions[:,:2] //= bases_scale
        composite.boxes.sizes[:,:2] *= phenotype_scale
        composite.boxes.sizes[:,:2] //= bases_scale

        # move boxes to be relative to this table
        composite.boxes.positions -= composite.boxes[index].position

        table = pandas.read_csv(input.table, index_col=0)
        #table['bbox_x1'] += composite.boxes[index].position[0]
        #table['bbox_y1'] += composite.boxes[index].position[1]
        #table['bbox_x2'] += composite.boxes[index].position[0]
        #table['bbox_y2'] += composite.boxes[index].position[1]
        #table.cells.bboxes += [[*composite.boxes[index].position, *composite.boxes[index].position]]
        #table['xpos'] += composite.boxes[index].position[0]
        #table['ypos'] += composite.boxes[index].position[1]
        #centroids = np.array([table['xpos'], table['ypos']]).T

        #boxes = constitch.BBoxList.from_table(table)

        neighbors = sklearn.neighbors.NearestNeighbors(n_neighbors=1).fit(composite.boxes.centers)
        distances, indices = neighbors.kneighbors(table.cells.centers)

        mask = indices == index
        debug (np.unique(indices, return_counts=True))
        debug ('mask ', mask.sum(), len(table.index))
        debug (composite.boxes[index].center)
        debug (table.cells.centers.mean(axis=0))

        max_cell_index = 0
        debug (table)

        for path in input.edge_tables:
            cur_table = pandas.read_csv(path, index_col=0)

            x, y = path.split('/tile')[1].split('y')[0].split('x')
            cur_index = int(x) * grid_size + int(y)

            cur_table['bbox_x1'] += composite.boxes[cur_index].position[0]
            cur_table['bbox_y1'] += composite.boxes[cur_index].position[1]
            cur_table['bbox_x2'] += composite.boxes[cur_index].position[0]
            cur_table['bbox_y2'] += composite.boxes[cur_index].position[1]

            debug (cur_table)
            overlapping_cells = table.cells.intersecting_cells(cur_table)
            debug (overlapping_cells)
            debug (list(overlapping_cells.index))
            #debug (overlapping_cells['area_ratio'])
            sums = {}
            for i,j in overlapping_cells.index:
                sums[i] = sums.get(i, 0) + overlapping_cells.cells[i,j].area()

            to_remove = [i for i, total in sums.items() if total > table.cells[i].area() * overlap_threshold]
            debug ('to_remove', len(to_remove))
            for i in to_remove:
                mask[table.index.get_loc(i)] = False
            """
            #cur_boxes = constitch.BBoxList.from_table(cur_table)
            max_cell_index = max(max_cell_index, max(cur_table.index))

            largest_cell = max(table.cells.sizes.max(), cur_table.cells.sizes.max())

            neighbors = sklearn.neighbors.NearestNeighbors(n_neighbors=5).fit(cur_boxes.centers)
            distances, indices = neighbors.radius_neighbors(boxes.centers, radius=largest_cell)
            debug ('found neighbors', distances.shape)

            for i in range(len(distances)):
                if not mask[i]: continue

                total_overlap = sum(boxes[i].intersection(cur_boxes[j]).area)
                mask[i] = total_overlap < table['area'].iloc[i] * overlap_threshold
            """

        debug ('mask ', mask.sum(), len(table.index))

        table = table[mask]
        table = table.reset_index(names='orig_index')
        table = table.set_index(pandas.RangeIndex(max_cell_index + 1, max_cell_index + 1 + len(table.index)))

        table.to_csv(output.table)





if config['segmentation'].get('match_masks', False):
    mask_pair = config['segmentation']['match_masks']
    if mask_pair is True:
        mask_pair = ['nuclei', 'cells']
    #print ('mask_pair', mask_pair)

    def grid_index_reference(wildcards):
        if '_grid' not in wildcards.path:
            return []

        groups = re.match('^(.*)_grid(\d+)/tile(\d+)x(\d+)y(.*)', wildcards.path).groups()
        path1, path2 = groups[0], groups[-1]
        grid_size, x, y = map(int, groups[1:-1])
        index = x * grid_size + y
        if index == 0:
            return []
        x, y = (index - 1) // grid_size, (index - 1) % grid_size
        newpath = '{}_grid{}/tile{:02}x{:02}y{}'.format(path1, grid_size, x, y, path2)

        paths = expand(segmentation_dir + newpath + '/{segmentation_type}{extra_params}.csv', segmentation_type=mask_pair, allow_missing=True)
        return paths

    rule match_cell_tables:
        """ Matches segmentations that should be labeled the same, typically cells and nuclei.
        Uses one of the segmentations as a base segmentation (normally nuclei as nuclei segmentation is more reliable),
        then the other segmentations are matched to the base segmentation. For each mask in the base segmentation
        the mask in the other segmentation with the greatest overlapping area is chosen. Any masks that don't have
        a mapping are removed
        """
        input:
            tables = expand(segmentation_dir + '{path}/{segmentation_type}{extra_params}_unmatched.csv', segmentation_type=mask_pair, allow_missing=True),
            grid_index_reference = grid_index_reference,
        output:
            tables = expand(segmentation_dir + '{path}/{segmentation_type}{extra_params}.csv', segmentation_type=mask_pair, allow_missing=True),
        wildcard_constraints:
            extra_params = '(|bases)(|expand\d+)',
        run:
            import numpy as np
            import pandas
            import constitch
            import starcall.cells

            base_table = pandas.read_csv(input.tables[0], index_col=0)
            masks_scale = base_table.cells.best_masks_scale
            all_tables = []

            overlap_threshold = 0.25

            mapping = np.full((len(base_table.index), len(input.tables) - 1), -1)

            for i in range(len(input.tables) - 1):
                table = pandas.read_csv(input.tables[i+1], index_col=0)
                all_tables.append(table)

                if table.cells.best_masks_scale != masks_scale:
                    table = table.copy()
                    table.cells.rescale_masks(masks_scale)

                intersecting = base_table.cells.intersecting_cells(table)
                cur_mapping = np.full(len(table.index), -1)
                ratios = np.zeros(len(table.index))

                debug (intersecting)
                intersecting = intersecting.sort_values('area', ascending=False)
                debug (intersecting)
                best_mapping = intersecting[~intersecting.index.to_frame().duplicated(0)]
                best_mapping = best_mapping[best_mapping['area_ratio']>=overlap_threshold]
                debug (best_mapping)
                reverse_mapping = intersecting[~intersecting.index.to_frame().duplicated(1)]
                reverse_mapping = reverse_mapping[reverse_mapping['area_ratio']>=overlap_threshold]
                debug (reverse_mapping)

                full_mapping = set(best_mapping.index) & set(reverse_mapping.index)
                for j, k in full_mapping:
                    mapping[base_table.index.get_loc(j),i] = k

                """
                for j, group in intersecting.groupby(level=0):
                    if len(group.index) == 0: continue

                    debug (group)
                    #group.to_csv('tmp.csv')
                    areas = np.array([group.cells[i].area() for i in group.index])
                    best = np.argmax(areas)
                    if group['area_ratio'].iloc[best] >= overlap_threshold:
                        mapping[base_table.index.get_loc(j),i] = group.index[best][1]

                largest_cell = max(table.cells.sizes.max(), base_table.cells.sizes.max())

                neighbors = sklearn.neighbors.NearestNeighbors(n_neighbors=5).fit(table.cells.centers)
                distances, indices = neighbors.radius_neighbors(base_table.cells.centers, radius=largest_cell)

                for j in range(len(indices)):
                    overlaps = [base_table.cells[j].intersection(table.cells[k]).area() for k in indices[j]]
                    best_index = np.argmax(overlaps)
                    if overlaps[best_index] / base_table.cells[j].area() >= overlap_threshold:
                        mapping[j,i] = table.index[indices[best_index]]
                """

            mask = np.all(mapping != -1, axis=1)
            base_table = base_table[mask]
            mapping = mapping[mask]

            # relabeling to be sequential

            max_indices = [0] * len(input.tables)
            # if in a grid, have to read previous tile to see max index
            if len(input.grid_index_reference):
                max_indices = [max(pandas.read_csv(path, index_col=0).index) for path in input.grid_index_reference]

            if 'orig_index' not in base_table.columns:
                base_table = base_table.reset_index(names='orig_index')
            base_table = base_table.set_index(pandas.RangeIndex(max_indices[0] + 1, max_indices[0] + 1 + len(base_table.index)))
            base_table.to_csv(output.tables[0])

            for i, table in enumerate(all_tables):
                table = table.loc[mapping[:,i],:]
                if 'orig_index' not in table.columns:
                    table = table.reset_index(names='orig_index')
                table = table.set_index(pandas.RangeIndex(max_indices[0] + 1, max_indices[i+1] + 1 + len(table.index)))
                table.to_csv(output.tables[i+1])

    ruleorder: match_cell_tables > split_grid_table


def get_segmentation_grid(wildcards):
    grid_size = int(wildcards.grid_size)
    numbers = ['{:02}'.format(i) for i in range(grid_size)]
    return expand(segmentation_dir + '{well}_grid{grid_size}/tile{x}x{y}y/{segmentation_type}.csv', x=numbers, y=numbers, allow_missing=True)

rule concat_cell_tables:
    """ Combine cell tables from grid tiles together into a single table. Because
    duplicate cells have been removed by drop_duplicate_cells, the tables can simply
    be concatenated together.
    """
    input:
        tables = get_segmentation_grid,
        composite = segmentation_dir + '{well}_grid{grid_size}/grid_composite.json',
    output:
        table = segmentation_dir + '{well}_grid{grid_size,\d+}/{segmentation_type}.csv',
    run:
        import pandas
        import numpy as np
        import constitch

        composite = constitch.load(input.composite)

        composite.boxes.positions[:,:2] *= phenotype_scale
        composite.boxes.positions[:,:2] //= bases_scale
        composite.boxes.sizes[:,:2] *= phenotype_scale
        composite.boxes.sizes[:,:2] //= bases_scale

        tables = []
        for i, path in enumerate(input.tables):
            table = pandas.read_csv(path, index_col=0)
            box = composite.boxes[i]

            table['bbox_x1'] += box.position[0]
            table['bbox_x2'] += box.position[0]
            table['bbox_y1'] += box.position[1]
            table['bbox_y2'] += box.position[1]
            table['orig_file_index'] = np.full(len(table.index), i)

            tables.append(table)

        table = pandas.concat(tables)
        table.to_csv(output.table)


### Merging then resplitting cell segmentation


def stitch_segmentation_section(image_paths, composite, section_box, table, phenotype=True):
    import constitch
    import numpy as np
    import tifffile
    import pandas

    if type(composite) == str:
        composite = constitch.load(composite)

    if type(table) == str:
        table = pandas.read_csv(table, index_col=0)

    if phenotype:
        # scaling from base images to phenotype
        composite.boxes.positions[:,:2] *= phenotype_scale
        composite.boxes.positions[:,:2] //= bases_scale
        composite.boxes.sizes[:,:2] *= phenotype_scale
        composite.boxes.sizes[:,:2] //= bases_scale

    dtype = [dtype for dtype in [np.uint16, np.uint32, np.uint64] if np.iinfo(dtype).max > len(table.index) + 1][0]
    touching_indices = [i for i, box in enumerate(composite.boxes) if box.collides(section_box)]

    #composite.images = [tifffile.memmap(path, mode='r') for path in image_paths]

    composite.images = []

    table['new_index'] = np.arange(1, len(table.index) + 1)
    debug (table)
    #second_mapping = {table.index[i]: i+1 for i in range(len(table.index))}
    #debug (second_mapping)
    #debug ('59667 mapped to', second_mapping.get(59667, 0))

    for i in range(len(composite.boxes)):
        if i not in touching_indices:
            composite.images.append(np.empty((1, 1), dtype))
            continue

        image = tifffile.imread(image_paths[i])
        cur_table = table[table['orig_file_index']==i]

        #mapping_arr = np.zeros(cur_table['orig_index'].max() + 1, dtype)
        mapping_arr = np.zeros(image.max() + 1, dtype)
        for orig_cell, new_cell in zip(cur_table['orig_index'], cur_table['new_index']):
            mapping_arr[orig_cell] = new_cell
        #mapping = mapping_table[mapping_table['table_index']==i]
        #for cell, newcell in zip(mapping['cell'], mapping['new_cell']):
            #mapping_arr[cell] = second_mapping.get(newcell, 0)

        #debug (mapping_arr.shape, image.max())
        #debug (max(second_mapping.keys()), max(second_mapping.values()))
        #debug (mapping['cell'].max(), mapping['new_cell'].max())
        image = mapping_arr[image]
        debug (i, 'is 18 in image', 18 in list(np.unique(image)))
        debug (i, 'is 17 in image', 17 in list(np.unique(image)))
        composite.images.append(image)

    debug (len(composite.boxes), len(composite.images), len(image_paths))
    #merger = constitch.EfficientNearestMerger()
    merger = constitch.MaxMerger()

    full_image = composite.stitch(merger=merger, indices=touching_indices, mins=section_box.point1, maxes=section_box.point2)
    del composite

    max_label, num_unique = full_image.max(), np.unique(full_image).shape[0]
    debug ('Max label', max_label, 'Num unique', num_unique)
    #assert max_label == num_unique - 1
    if max_label != num_unique - 1:
        debug ('BIG PROBLEM cellprofiler will not like this')
        debug (set(range(max_label + 1)) - set(np.unique(full_image)))
        assert max_label == num_unique - 1

        """
        missing_index = next(iter(set(range(max_label + 1)) - set(np.unique(full_image))))
        nuclei = tifffile.imread('segmentation/tmpwell3_grid20/tile12x14y/nuclei_mask.tif')

        for index in [missing_index, 1016, 1891]:
            debug ('index', index)
            orig_index = table.index[index-1]
            xpos, ypos = int(table['bbox_x1'][orig_index]), int(table['bbox_y1'][orig_index])
            radius = 500
            section_x1, section_y1 = max(0, xpos - radius), max(0, ypos - radius)
            section = full_image[section_x1:xpos+radius,section_y1:ypos+radius]
            section2 = nuclei[section_x1:xpos+radius,section_y1:ypos+radius]
            section, section2 = section.astype(int), section2.astype(int)
            section[section==0] = -1000
            section2[section2==0] = -1000
            debug (xpos, ypos, section_x1, section_y1, section.shape)

            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(nrows=2, figsize=(5, 8))
            axes[0].imshow(section)
            axes[1].imshow(section2)
            x1, y1, x2, y2 = table['bbox_x1'][orig_index], table['bbox_y1'][orig_index], table['bbox_x2'][orig_index], table['bbox_y2'][orig_index]
            x1, y1, x2, y2 = x1 - section_x1, y1 - section_y1, x2 - section_x1, y2 - section_y1
            debug (x1, y1, x2, y2)
            debug (np.unique(section[x1:x2,y1:y2]))
            debug (np.unique(section2[x1:x2,y1:y2]))
            axes[0].plot([y1, y2, y2, y1, y1], [x1, x1, x2, x2, x1])
            axes[1].plot([y1, y2, y2, y1, y1], [x1, x1, x2, x2, x1])
            fig.savefig('tmp_missing_cell{}.png'.format(index))
        """

    return full_image


unmatched = '_unmatched' if config['segmentation'].get('match_masks', False) else ''

rule relabel_segmentation:
    """ Uses a filtered and matched cells table to relabel a segmentation mask to be sequential.
    """
    input:
        image = segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/{segmentation_type}_mask{downscaled}' + unmatched + '{grid}.tif',
        table = segmentation_dir + '{path_nogrid}{grid}{path_nogrid2}/{segmentation_type}.csv',
    output:
        image = '{output_dir}{path_nogrid}{grid}{path_nogrid2}/{segmentation_type}_mask{downscaled,|_downscaled}.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 5 + 10000,
    run:
        import tifffile
        import pandas
        import numpy as np

        image = tifffile.imread(input.image)
        table = pandas.read_csv(input.table, index_col=0)

        dtype = [dtype for dtype in [np.uint16, np.uint32, np.uint64] if np.iinfo(dtype).max > len(table.index) + 1][0]
        mapping_arr = np.zeros(image.max() + 1, dtype)
        for i, orig_cell in enumerate(table['orig_index']):
            mapping_arr[orig_cell] = i + 1

        image = mapping_arr[image]
        tifffile.imwrite(output.image, image)

ruleorder: relabel_segmentation > stitch_tile_segmentation



segmentation_grid_size = config.get('segmentation_grid_size', 1)

def find_othertable(wildcards):
    match_pair = config['segmentation'].get('match_masks', False)
    if match_pair:
        for segtype in match_pair[1:]:
            if wildcards.segmentation_type.count(segtype) != 0:
                newtype = wildcards.segmentation_type.replace(segtype, match_pair[0])
                return ['{output_dir}{well}_grid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/' + newtype + '.csv']
    #if config['segmentation'].get('match_masks', False) and wildcards.segmentation_type.count('cells') != 0:
        #newtype = wildcards.segmentation_type.replace('cells', 'nuclei')
        #return ['{output_dir}{well}_grid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/' + newtype + '.csv']
    return []


rule split_grid_table:
    """ Creates a cell table for a grid tile from a full well cell table, including only cells that
    are closest to the given tile. This way each cell has a single tile it is included in, and there will
    be no duplicate cells when the tables are rejoined together.

    The only complication to this is when splitting tables that have been matched to each other, eg cell and nuclei
    tables. In this case it is not guaranteed that the cell and nucleus with the same label will be mapped to the
    same tile. To resolve this, only one of the matched tables is split, eg the nuclei table. Then the cell table
    is split using the same indices as the nuclei table. This is accomplished with input.othertable, which is empty
    unless the table is a type of segmentation that is matched, in which case it will be the matched table.
    """
    input:
        #table = segmentation_dir + '{well}_grid/{segmentation_type}.csv',
        table = segmentation_dir + '{well}_grid' + str(segmentation_grid_size) + '/{segmentation_type}.csv',
        composite = segmentation_dir + '{well}_grid{grid_size}/grid_composite.json',
        othertable = find_othertable,
    output:
        table = '{output_dir}{well}_grid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/{segmentation_type}.csv'
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 2
    run:
        import pandas
        import numpy as np
        import constitch
        import sklearn.neighbors
        import starcall.cells

        table = pandas.read_csv(input.table, index_col=0)
        composite = constitch.load(input.composite)
        grid_size, x, y = int(wildcards.grid_size), int(wildcards.x), int(wildcards.y)
        index = x * grid_size + y

        composite.boxes.positions[:,:2] *= phenotype_scale
        composite.boxes.positions[:,:2] //= bases_scale
        composite.boxes.sizes[:,:2] *= phenotype_scale
        composite.boxes.sizes[:,:2] //= bases_scale
        box = composite.boxes[index]
        debug (table)

        if len(input.othertable) == 0:
            neighbors = sklearn.neighbors.NearestNeighbors(n_neighbors=1).fit(composite.boxes.centers)
            distances, indices = neighbors.kneighbors(table.cells.centers)
            table = table.loc[indices==index]
        else:
            othertable = pandas.read_csv(input.othertable[0], index_col=0)
            table = table.loc[othertable.index]

        table = table.copy()
        table['bbox_x1'] -= box.position[0]
        table['bbox_x2'] -= box.position[0]
        table['bbox_y1'] -= box.position[1]
        table['bbox_y2'] -= box.position[1]
        debug (table)

        table.to_csv(output.table)



def get_grid_filenames(wildcards):
    grid_size = int(wildcards.grid_size)
    numbers = ['{:02}'.format(i) for i in range(segmentation_grid_size)]
    grid = '_grid' + str(segmentation_grid_size)
    unmatched = '_unmatched' if config['segmentation'].get('match_masks', False) else ''
    return expand(segmentation_dir + '{well}' + grid
                + '/tile{x}x{y}y/{segmentation_type}_mask{downscaled}' + unmatched
                + grid + '.tif', x=numbers, y=numbers, allow_missing=True)

def get_cells_mapping2(wildcards):
    if config['segmentation'].get('match_masks', False):
        return (segmentation_dir + '{well}_grid' + str(segmentation_grid_size) + '/'
                + wildcards.segmentation_type.replace('cells', 'nuclei')
                + '_mappings.csv')
    return segmentation_dir + '{well}_grid{grid_size,\d+}/{segmentation_type}_mappings.csv'

rule stitch_tile_segmentation:
    """ Stitches segmentation masks from one grid of tiles to another grid of tiles. The main difficulty
    in this is making sure the masks are relabeled properly, ensuring the resulting mask is labeled sequentially.
    The cell table corresponding to the output segmentation is used as a key, with masks labeled according to
    the order in which they appear in the table (eg the first row will be labeled as 1, the second as 2 ...)

    Because the cell table contains a column 'orig_index' it is easy to know the label of each cell in the input masks,
    and they are relabeled before stitching.
    """
    input:
        images = get_grid_filenames,
        composite = segmentation_dir + '{well}_grid' + str(segmentation_grid_size) + '/grid_composite.json',
        #mappings = segmentation_dir + '{well}_cellgrid' + str(segmentation_grid_size) + '/{segmentation_type}_mappings.csv',
        #mappings = get_cells_mapping2,
        composite2 = stitching_dir + '{well}_grid{grid_size}/grid_composite.json',
        table = '{output_dir}{well}_grid{grid_size}/tile{x}x{y}y/{segmentation_type}.csv',
    output:
        image = temp('{output_dir}{well}_grid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/{segmentation_type}_mask{downscaled,|_downscaled}.tif'),
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 2
    run:
        import tifffile
        import constitch

        composite = constitch.load(input.composite2)

        if wildcards.downscaled == '':
            # scaling from base images to phenotype
            composite.boxes.positions[:,:2] *= phenotype_scale
            composite.boxes.positions[:,:2] //= bases_scale
            composite.boxes.sizes[:,:2] *= phenotype_scale
            composite.boxes.sizes[:,:2] //= bases_scale

        tile_index = int(wildcards.x) * int(wildcards.grid_size) + int(wildcards.y)

        tifffile.imwrite(output.image, stitch_segmentation_section(input.images,
                input.composite, composite.boxes[tile_index], input.table, phenotype=wildcards.downscaled == ''))

"""
def get_cells_mapping3(wildcards):
    if config['segmentation'].get('match_masks', False):
        return (segmentation_dir + '{well}_grid{grid_size}/'
                + wildcards.segmentation_type.replace('cells', 'nuclei')
                + '_mappings.csv')
    return segmentation_dir + '{well}_cellgrid{grid_size,\d+}/{segmentation_type}_mappings.csv'

rule stitch_well_segmentation:
    input:
        images = get_grid_filenames,
        composite = segmentation_dir + '{well}_grid{grid_size}/grid_composite.json',
        mappings = get_cells_mapping3,
        table = '{output_dir}{well}_grid{grid_size}/{segmentation_type}.csv',
    output:
        image = '{output_dir}{well}_grid{grid_size,\d+}/{segmentation_type}_mask.tif',
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 2
    run:
        import tifffile
        import constitch

        composite = constitch.load(input.composite)
        # scaling from base images to phenotype
        composite.boxes.positions[:,:2] *= phenotype_scale
        composite.boxes.positions[:,:2] //= bases_scale
        composite.boxes.sizes[:,:2] *= phenotype_scale
        composite.boxes.sizes[:,:2] //= bases_scale

        mins, maxes = composite.boxes.points1.min(axis=0)[:2], composite.boxes.points2.max(axis=0)[:2]

        tifffile.imwrite(output.image, stitch_segmentation_section(input.images,
                input.composite, input.mappings, constitch.BBox(point1=mins, point2=maxes), input.table))


rule stitch_tile_from_well_segmentation:
    input:
        image = segmentation_dir + '{well}_grid/{segmentation_type}_mask_unmerged.tif',
        composite = segmentation_dir + '{well}_grid{grid_size}/grid_composite.json',
        table = '{output_dir}{well}_grid{grid_size}/tile{x}x{y}y/{segmentation_type}.csv',
        full_table = segmentation_dir + '{well}_grid/{segmentation_type}_unmerged.csv',
    output:
        image = '{output_dir}{well}_grid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/{segmentation_type}_mask.tif',
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 2
    run:
        import tifffile
        import constitch
        import numpy as np
        import pandas

        image = tifffile.memmap(input.image, mode='r')
        input_composite = constitch.Composite([image])
        table = pandas.read_csv(input.full_table, index_col=0)
        fake_mapping = pandas.DataFrame(dict(
                table_index=np.full(len(table.index), 0),
                cell=table.index, new_cell=table.index))

        composite = constitch.load(input.composite)
        # scaling from base images to phenotype
        composite.boxes.positions[:,:2] *= phenotype_scale
        composite.boxes.positions[:,:2] //= bases_scale
        composite.boxes.sizes[:,:2] *= phenotype_scale
        composite.boxes.sizes[:,:2] //= bases_scale
        tile_index = int(wildcards.x) * int(wildcards.grid_size) + int(wildcards.y)

        tifffile.imwrite(output.image, stitch_segmentation_section([input.image],
                input_composite, fake_mapping, composite.boxes[tile_index], input.table))
"""




def get_grid_size_file(wildcards):
    grid_size = config.get('segmentation_{}_grid_size'.format(wildcards.segmentation_type), segmentation_grid_size)
    if grid_size == 1:
        return segmentation_dir + '{well}/{segmentation_type}{filetype}'
    return segmentation_dir + '{well}_grid' + str(grid_size) + '/{segmentation_type}{filetype}',

rule link_merged_grid:
    input:
        get_grid_size_file,
    output:
        segmentation_dir + '{well}_grid/{segmentation_type}{filetype,.csv|_mask.tif}',
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"

ruleorder: link_merged_grid > segment_cells
ruleorder: link_merged_grid > segment_nuclei
ruleorder: link_merged_grid > tabulate_cells









#### QC

rule make_cell_overlay:
    """ Overlays the cell segmentation boundaries onto the phenotyping images used to
    create the segmentation. Useful to make sure segmentation is working well and to
    test different parameters, such as the diameter provided to cellpose.
    """
    input:
        image = (segmentation_dir + '{path}/corrected_pt.tif'
                if config['segmentation']['use_corrected'] else
                segmentation_dir + '{path}/raw_pt.tif'),
        cells = segmentation_dir + '{path}/{segmentation_type}_mask{params}.tif',
    output:
        qc_dir + '{path}/{segmentation_type}_overlay{params}.tif',
        qc_dir + '{path}/{segmentation_type}_overlay{params}.png',
    wildcard_constraints:
        params = params_regex('diameter', 'nuclearchannel', 'cytochannel'),
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 15 + 10000
    run:
        import numpy as np
        import tifffile
        import starcall.utils

        cells = tifffile.imread(input.cells)
        cell_borders = np.roll(cells, (1,1), axis=(0,1)) != cells
        debug (cell_borders.min(), cell_borders.max())
        cell_borders = cell_borders | (np.roll(cells, (-1,1), axis=(0,1)) != cells)
        debug (cell_borders.min(), cell_borders.max())
        cell_mask = cells != 0
        del cells

        tmp_image = tifffile.memmap(input.image, mode='r')
        debug(tmp_image.shape)
        image = np.zeros((tmp_image.shape[1] + 1, *tmp_image.shape[2:]), tmp_image.dtype)
        del tmp_image
        out_image = image[:-1]
        tifffile.imread(input.image, out=out_image.reshape(1, *out_image.shape))
        image[-1] = cell_borders * np.iinfo(image.dtype).max
        np.maximum(image[-1], cell_mask * image.dtype.type(np.iinfo(image.dtype).max // 3), out=image[-1])
        debug (image[-1].min(), image[-1].max())
        tifffile.imwrite(output[0], image)
        rgbimage = starcall.utils.to_rgb8(image)
        tifffile.imwrite(output[1], rgbimage)

