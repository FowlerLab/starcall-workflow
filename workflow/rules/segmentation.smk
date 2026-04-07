import os
import glob
import re


def get_segmentation_pt(wildcards):
    path = wildcards.path_nogrid.replace('_cellgrid', '_grid')
    if config['segmentation']['use_corrected']:
        return segmentation_dir + path + '/corrected_pt.tif'
    else:
        return segmentation_dir + path + '/raw_pt.tif'

rule segment_nuclei:
    """ Uses Stardist to segment the nuclei of cells in the phenotyping
    images. Takes one phenotyping channel as a nuclear channel. Outputs
    an image with integer masks for each cell.

    Params:
        nuclearchannel: The phenotyping channel to perform segmentation on, should be
            an integer index or one of the phenotyping channels specified in config.yaml
    """
    input:
        get_segmentation_pt,
    output:
        segmentation_dir + '{path_nogrid}/nuclei{nuclearchannel}_mask_unmatched.tif',
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
            nuclei = skimage.segmentation.clear_border(nuclei)
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
        get_segmentation_pt,
    output:
        segmentation_dir + '{path_nogrid}/cells{diameter}{nuclearchannel}{cytochannel}_mask_unmatched.tif',
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
            cells = skimage.segmentation.clear_border(cells)

            debug ('Found', cells.max(), 'cells')

            tifffile.imwrite(output[0], cells)#, compression='deflate')


rule expand_segmentation:
    input:
        cells = segmentation_dir + '{path_nogrid}/{segmentation_type}_mask.tif',
    output:
        cells = segmentation_dir + '{path_nogrid}/{segmentation_type}expanded{size,\d+}_mask.tif',
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


def get_segmentation_bases(wildcards):
    path = wildcards.path_nogrid.replace('_cellgrid', '_grid')
    return segmentation_dir + path + '/raw.tif'

rule segment_cells_bases:
    input:
        get_segmentation_bases,
    output:
        segmentation_dir + '{path_nogrid}/cellsbases{diameter}{nuclearchannel}_mask_unmatched.tif',
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

            debug(f'found {cells.max()} cells ')
            tifffile.imwrite(output[1], cells)#, compression='deflate')


rule segment_nuclei_bases:
    input:
        get_segmentation_bases,
    output:
        segmentation_dir + '{path_nogrid}/nucleibases{nuclearchannel}_mask_unmatched.tif',
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
            debug(f'found {nuclei.max()} nuclei ')

            tifffile.imwrite(output[0], nuclei)#, compression='deflate')


rule downscale_segmentation:
    """ Rescales cell segmentation to match the scale of the sequencing images. The
    scale ratio is specified in config.yaml, with phenotype_scale and sequencing_scale.
    """
    input:
        segmentation_dir + '{path_nogrid}/{segmentation_type}_mask.tif',
    output:
        segmentation_dir + '{path_nogrid}/{segmentation_type}_mask_downscaled.tif',
    run:
        import tifffile
        import skimage.transform

        mask = tifffile.imread(input[0])

        if wildcards.segmentation_type.count('bases') != 0:
            tifffile.imwrite(output[0], mask)
        else:
            tifffile.imwrite(output[0], skimage.transform.rescale(mask, bases_scale/phenotype_scale, order=0))


rule make_cell_overlay:
    """ Overlays the cell segmentation boundaries onto the phenotyping images used to
    create the segmentation. Useful to make sure segmentation is working well and to
    test different parameters, such as the diameter provided to cellpose.
    """
    input:
        image = get_segmentation_pt,
        cells = segmentation_dir + '{path_nogrid}/{segmentation_type}_mask{params}.tif',
    output:
        qc_dir + '{path_nogrid,.*}/{segmentation_type}_overlay{params}.tif',
        qc_dir + '{path_nogrid,.*}/{segmentation_type}_overlay{params}.png',
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


if config['segmentation'].get('match_masks', False):
    rule match_masks:
        """ Links cells and nuclei segmented previously, by matching overlapping
        cells/nuclei and choosing the pair with the highest overlap. Cells with no
        nuclei and nuclei with no cells are discarded.
        """
        input:
            cells = segmentation_dir + '{path_nogrid}/cells_mask_unmatched.tif',
            nuclei = segmentation_dir + '{path_nogrid}/nuclei_mask_unmatched.tif',
        output:
            cells = segmentation_dir + '{path_nogrid}/cells_mask.tif',
            nuclei = segmentation_dir + '{path_nogrid}/nuclei_mask.tif',
        resources:
            mem_mb = lambda wildcards, input: input.size_mb * 5 + 5000
        run:
            import collections
            import tifffile
            import numpy as np
            import skimage.measure
            import starcall.segmentation

            cells = tifffile.imread(input.cells)
            nuclei = tifffile.imread(input.nuclei)

            cells, nuclei = starcall.segmentation.match_segmentations(cells, nuclei)

            tifffile.imwrite(output.cells, cells)
            tifffile.imwrite(output.nuclei, nuclei)
else:
    rule match_masks:
        input:
            segmentation_dir + '{path_nogrid}/{segmentation_type}_mask_unmatched.tif',
        output:
            segmentation_dir + '{path_nogrid}/{segmentation_type}_mask.tif',
        localrule: True
        shell:
            "cp -l {input[0]} {output[0]}"


rule tabulate_cells:
    """ Simple information is recorded about the segmented cells, such as position, bbox.
    """
    input:
        cells = segmentation_dir + '{path_nogrid}/{segmentation_type}_mask.tif',
    output:
        table = segmentation_dir + '{path_nogrid}/{segmentation_type}.csv',
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 1.5
    run:
        import numpy as np
        import tifffile
        import skimage.measure
        import starcall.sequencing
        import pandas

        cells = tifffile.imread(input.cells)

        table = {}

        props = skimage.measure.regionprops(cells, cache=False)
        props = {prop.label: prop for prop in props}

        index = sorted(list(props.keys()))

        bboxes = np.array([props[cell].bbox for cell in index])
        centroids = np.array([props[cell].centroid for cell in index])

        table['xpos'] = centroids[:,0]
        table['ypos'] = centroids[:,1]
        table['bbox_x1'] = bboxes[:,0]
        table['bbox_y1'] = bboxes[:,1]
        table['bbox_x2'] = bboxes[:,2]
        table['bbox_y2'] = bboxes[:,3]

        table = pandas.DataFrame(table, index=index)
        table.to_csv(output.table)



### Merging then resplitting cell segmentation

def get_grid_filenames(wildcards):
    grid_size = int(wildcards.grid_size)
    numbers = ['{:02}'.format(i) for i in range(grid_size)]
    return expand(segmentation_dir + '{well}_cellgrid{grid_size}/tile{x}x{y}y/{segmentation_type}_mask_unmatched.tif', x=numbers, y=numbers, allow_missing=True)

'''
rule merge_grid_segmentation:
    """ If cell segmentation was run on a grid of tiles, the segmentations must be merged together.
    For each cell mask in an overlapping region between tiles, a possible match is found based on overlap.
    If one mask overlaps the other to a high degree (>75%) the masks are merged into one. This process
    is somewhat error prone, so it is encouraged to keep the tiles used for cell segmentation large, so
    not many cells have to be merged.
    """
    input:
        images = get_grid_filenames,
        composite = stitching_dir + '{well}_grid{grid_size}/grid_composite.json',
    output:
        image = segmentation_dir + '{well}_cellgrid{grid_size,\d+}/{segmentation_type}_mask_unmatched3.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 10 + 25000
    run:
        import numpy as np
        import tifffile
        import constitch
        import pandas

        composite = constitch.load(input.composite)

        num_cells = 0
        for i,path in enumerate(input.images):
            composite.images[i] = tifffile.imread(path)
            num_cells += composite.images[i].max()

        dtype = [dtype for dtype in [np.uint16, np.uint32, np.uint64] if np.iinfo(dtype).max > num_cells + 1][0]

        # scaling from base images to phenotype
        composite.boxes.positions[:,:2] *= phenotype_scale
        composite.boxes.positions[:,:2] //= bases_scale
        composite.boxes.sizes[:,:2] *= phenotype_scale
        composite.boxes.sizes[:,:2] //= bases_scale
        merger = constitch.EfficientMaskMerger(dtype=dtype)

        full_image = composite.stitch(merger=merger)
        del composite
        
        max_label, num_unique = full_image.max(), np.unique(full_image).shape[0]
        debug ('Max label', max_label, 'Num unique', num_unique)
        assert max_label == num_unique - 1
        tifffile.imwrite(output.image, full_image)
'''

def get_segmentation_grid(wildcards):
    grid_size = int(wildcards.grid_size)
    numbers = ['{:02}'.format(i) for i in range(grid_size)]
    return expand(segmentation_dir + '{well}_cellgrid{grid_size}/tile{x}x{y}y/{segmentation_type}.csv', x=numbers, y=numbers, allow_missing=True)

rule merge_segmentation_tables:
    input:
        tables = get_segmentation_grid,
        composite = stitching_dir + '{well}_grid{grid_size}/grid_composite.json',
    output:
        #table = segmentation_dir + '{well}_cellgrid{grid_size,\d+}/{segmentation_type}.csv',
        mappings = segmentation_dir + '{well}_cellgrid{grid_size,\d+}/{segmentation_type}_mappings.csv',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 10 + 10000
    run:
        import pandas
        import numpy as np
        import constitch
        import sklearn.neighbors
        import matplotlib.pyplot as plt

        overlap_threshold = 0.5

        composite = constitch.load(input.composite)

        # scaling from base images to phenotype
        composite.boxes.positions[:,:2] *= phenotype_scale
        composite.boxes.positions[:,:2] //= bases_scale
        composite.boxes.sizes[:,:2] *= phenotype_scale
        composite.boxes.sizes[:,:2] //= bases_scale

        tables = []
        boxes = []
        table_indices = []
        cell_indices = []
        for i,path in enumerate(input.tables):
            table = pandas.read_csv(path, index_col=0)
            tables.append(table)

            boxes.extend(constitch.BBox(
                    point1=[table['bbox_x1'][j] + composite.boxes[i].point1[0],
                        table['bbox_y1'][j] + composite.boxes[i].point1[1]],
                    point2=[table['bbox_x2'][j] + composite.boxes[i].point1[0],
                        table['bbox_y2'][j] + composite.boxes[i].point1[1]])
                        for j in table.index)

            for j in table.index:
                box = constitch.BBox(
                        point1=[table['bbox_x1'][j] + composite.boxes[i].point1[0],
                            table['bbox_y1'][j] + composite.boxes[i].point1[1]],
                        point2=[table['bbox_x2'][j] + composite.boxes[i].point1[0],
                            table['bbox_y2'][j] + composite.boxes[i].point1[1]])
                if box.size.max() > 1000:
                    debug (box.point1, box.point2)
                    debug (
                        'point1', [table['bbox_x1'][j] + composite.boxes[i].point1[0],
                            table['bbox_y1'][j] + composite.boxes[i].point1[1]],
                        'point2', [table['bbox_x2'][j] + composite.boxes[i].point1[0],
                            table['bbox_y2'][j] + composite.boxes[i].point1[1]])
                    debug ('tile box', composite.boxes[i].point1, composite.boxes[i].point2)
                    debug (table['bbox_x1'][j], table['bbox_y1'][j], table['bbox_x2'][j], table['bbox_y2'][j])
                    debug (i, j)
                    skjdflskjdflsk

            table_indices.extend([i] * len(table.index))
            cell_indices.extend(table.index)


        largest_cell = max(np.linalg.norm(box.size) for box in boxes)
        centers = np.array([box.center for box in boxes])
        debug ('centers', centers.shape, largest_cell)

        neighbors = sklearn.neighbors.NearestNeighbors(n_neighbors=5).fit(centers)
        distances, indices = neighbors.radius_neighbors(centers, radius=largest_cell)
        debug ('found neighbors', distances.shape)

        fig, axes = plt.subplots(figsize=(50,50))
        fig2, axes2 = plt.subplots(figsize=(50,50))
        fig3, axes3 = plt.subplots(figsize=(50,50))
        for i, ibox in enumerate(composite.boxes):
            ibox.plot(axes)
            ibox.plot(axes2)
            ibox.plot(axes3, color='C{}'.format(i%10))

        matches = {}
        matches_reverse = {}
        to_remove = set()
        for i in range(len(boxes)):
            should_be_matched = sum(box.intersection(boxes[i]).area() > 0 for box in composite.boxes) > 1

            dists_to_others = []
            for k, other_image_box in enumerate(composite.boxes):
                if k == table_indices[i]: continue
                #dist_to_edge = np.minimum(boxes[i].point2 - other_image_box.point1, other_image_box.point2 - boxes[i].point1).min()
                dist_to_edge = np.minimum(boxes[i].point1 - other_image_box.point1, other_image_box.point2 - boxes[i].point2).min()
                dists_to_others.append(dist_to_edge)

            image_box = composite.boxes[table_indices[i]]
            dist_to_edge = np.minimum(boxes[i].point2 - image_box.point1, image_box.point2 - boxes[i].point1).min()
            #dist_to_edge = np.minimum(box.point1 - image_box.point1, image_box.point2 - box.point2).min()
            #print (dists_to_others, dist_to_edge, file=sys.stderr)
            if any((dist > dist_to_edge) for dist in dists_to_others):
                #boxes[i].plot(axes3)
                debug ('too farr', boxes[i].point1, boxes[i].point2)
                to_remove.add(i)
                continue

            for j, dist in zip(indices[i], distances[i]):
                if table_indices[i] == table_indices[j]:
                    continue

                area = boxes[i].intersection(boxes[j]).area()
                if area / boxes[i].area() >= overlap_threshold and area / boxes[j].area() >= overlap_threshold:
                    pair = [i, j]
                    if pair[0] > pair[1]: pair = [pair[1], pair[0]]

                    while pair[0] in matches:
                        pair[0] = matches[pair[0]]

                    boxes[pair[0]].plot(axes, color='C{}'.format(i % 10))
                    boxes[pair[1]].plot(axes, color='C{}'.format(i % 10))

                    matches[pair[1]] = pair[0]
                    matches_reverse[pair[0]] = pair[1]
                    break
            else:
                if should_be_matched:
                    boxes[i].plot(axes2)
                    boxes[i].plot(axes, color='grey')
                #else:
                    #boxes[i].plot(axes2, color='grey')

        fig.savefig('tmp_boxes_matched_50p.png')
        fig2.savefig('tmp_boxes_matched2_50p.png')

        mapping = []
        label = 1
        for i in range(len(boxes)):
            if i in matches:
                mapping.append(mapping[matches[i]])
            elif i in to_remove and i not in matches_reverse:
                boxes[i].plot(axes3, color='C{}'.format(table_indices[i]%10))
                mapping.append(0)
            else:
                mapping.append(label)
                label += 1

        fig3.savefig('tmp_boxes_matched3_50p.png')

        mapping_table = pandas.DataFrame(dict(table_index=table_indices, cell=cell_indices, new_cell=mapping))
        mapping_table.to_csv(output.mappings)

        """
        full_table = pandas.concat(tables, ignore_index=True)
        full_table['cell'] = mapping
        full_table = full_table.drop_duplicates(subset='cell', keep='first')
        full_table = full_table.set_index('cell')
        #full_table = full_table.set_index(mapping)
        full_table.to_csv(output.table)

        dtype = [dtype for dtype in [np.uint16, np.uint32, np.uint64] if np.iinfo(dtype).max > len(boxes) + 1][0]

        # scaling from base images to phenotype
        composite.boxes.positions[:,:2] *= phenotype_scale
        composite.boxes.positions[:,:2] //= bases_scale
        composite.boxes.sizes[:,:2] *= phenotype_scale
        composite.boxes.sizes[:,:2] //= bases_scale

        merger = constitch.EfficientMaskMerger(dtype=dtype)

        mapping_table = dict(table_index=[], cell=[], new_cell=[])
        #filtered_tables = []

        for i, table in enumerate(tables):
            boxes = [constitch.BBox(
                    point1=[table['bbox_x1'][i] + composite.boxes[i].point1[0],
                        table['bbox_x2'][i] + composite.boxes[i].point2[0]],
                    point2=[table['bbox_y1'][i] + composite.boxes[i].point1[1],
                        table['bbox_y2'][i] + composite.boxes[i].point2[1]])
                        for i in table.index]
            mapping = merger.find_mapping(composite.boxes[i], boxes)

            mapping_table['table_index'].extend(i for j in range(len(table.index)))
            mapping_table['cell'].extend(table.index)
            mapping_table['new_cell'].extend(mapping)

        mapping_table = pandas.DataFrame(mapping_table).set_index('table_index')
        mapping_table.to_csv(output.mappings)
        """

def get_cells_mapping(wildcards):
    if config['segmentation'].get('match_masks', False):
        return (segmentation_dir + '{well}_cellgrid{grid_size,\d+}/'
                + wildcards.segmentation_type.replace('cells', 'nuclei')
                + '_mappings.csv')
    return segmentation_dir + '{well}_cellgrid{grid_size,\d+}/{segmentation_type}_mappings.csv'

rule merge_tables_mapping:
    input:
        tables = get_segmentation_grid,
        mappings = get_cells_mapping,
        composite = stitching_dir + '{well}_grid{grid_size}/grid_composite.json',
        #mappings = lambda wildcards: (segmentation_dir + '{well}_cellgrid{grid_size,\d+}/'
                #+ (wildcards.segmentation_type.replace('cells', 'nuclei')
                    #if config['segmentation'].get('match_masks', False) else wildcards.segmentation_type)
                #+ '_mappings.csv')
    output:
        table = segmentation_dir + '{well}_cellgrid{grid_size,\d+}/{segmentation_type}.csv',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 10 + 10000
    run:
        import pandas
        import numpy as np
        import constitch

        composite = constitch.load(input.composite)

        # scaling from base images to phenotype
        composite.boxes.positions[:,:2] *= phenotype_scale
        composite.boxes.positions[:,:2] //= bases_scale
        composite.boxes.sizes[:,:2] *= phenotype_scale
        composite.boxes.sizes[:,:2] //= bases_scale

        tables = []
        for i,path in enumerate(input.tables):
            table = pandas.read_csv(path, index_col=0)
            table['bbox_x1'] += composite.boxes[i].point1[0]
            table['bbox_y1'] += composite.boxes[i].point1[1]
            table['bbox_x2'] += composite.boxes[i].point1[0]
            table['bbox_y2'] += composite.boxes[i].point1[1]
            tables.append(table)

        mapping_table = pandas.read_csv(input.mappings)

        full_table = pandas.concat(tables, ignore_index=True)
        full_table['cell'] = mapping_table['new_cell']
        full_table = full_table.drop_duplicates(subset='cell', keep='first')
        full_table = full_table.set_index('cell')
        full_table = full_table.sort_index()
        full_table.to_csv(output.table)
        

ruleorder: merge_segmentation_tables > tabulate_cells

def stitch_segmentation_section(image_paths, composite_path, mappings_path, section_box, section_table_path):
    import constitch
    import numpy as np
    import tifffile
    import pandas

    composite = constitch.load(composite_path)

    # scaling from base images to phenotype
    composite.boxes.positions[:,:2] *= phenotype_scale
    composite.boxes.positions[:,:2] //= bases_scale
    composite.boxes.sizes[:,:2] *= phenotype_scale
    composite.boxes.sizes[:,:2] //= bases_scale

    table = pandas.read_csv(section_table_path, index_col=0)

    dtype = [dtype for dtype in [np.uint16, np.uint32, np.uint64] if np.iinfo(dtype).max > len(table.index) + 1][0]
    touching_indices = [i for i, box in enumerate(composite.boxes) if box.collides(section_box)]

    #composite.images = [tifffile.memmap(path, mode='r') for path in image_paths]

    composite.images = []
    mapping_table = pandas.read_csv(mappings_path)

    second_mapping = {table.index[i]: i+1 for i in range(len(table.index))}

    for i in range(len(composite.boxes)):
        if i not in touching_indices:
            composite.images.append(np.empty((1, 1), dtype))
            continue

        image = tifffile.imread(image_paths[i])
        mapping = mapping_table[mapping_table['table_index']==i]
        mapping_arr = np.zeros(mapping['cell'].max() + 1, dtype)
        for cell, newcell in zip(mapping['cell'], mapping['new_cell']):
            mapping_arr[cell] = second_mapping.get(newcell, 0)

        debug (mapping_arr.shape, image.max())
        debug (max(second_mapping.keys()), max(second_mapping.values()))
        debug (mapping['cell'].max(), mapping['new_cell'].max())
        image = mapping_arr[image]
        composite.images.append(image)

    debug (len(composite.boxes), len(composite.images), len(image_paths))
    merger = constitch.EfficientNearestMerger()

    full_image = composite.stitch(merger=merger, indices=touching_indices, mins=section_box.point1, maxes=section_box.point2)
    del composite
    
    max_label, num_unique = full_image.max(), np.unique(full_image).shape[0]
    debug ('Max label', max_label, 'Num unique', num_unique)
    #assert max_label == num_unique - 1

    return full_image






segmentation_grid_size = config.get('segmentation_grid_size', 1)

def get_grid_size_file(wildcards):
    grid_size = config.get('segmentation_{}_grid_size'.format(wildcards.segmentation_type), segmentation_grid_size)
    if grid_size == 1:
        return segmentation_dir + '{well}/{segmentation_type}{filetype}'
    return segmentation_dir + '{well}_cellgrid' + str(grid_size) + '/{segmentation_type}{filetype}',

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
ruleorder: link_merged_grid > match_masks
ruleorder: link_merged_grid > tabulate_cells
ruleorder: merge_tables_mapping > tabulate_cells


rule split_grid_table:
    """ Splits the cell info table into a tile in a grid, for use in later steps such as sequencing
    or phenotyping. Cells are matched to the tile closest to their centroid, ensuring no cells
    are included in two tiles. The overlap between tiles, specified in config.yaml, should be
    large enough that nearly all cells are contained in at least one tile. If enough cells are not
    contained in any tile (>20%) this step will fail and the overlap should be increased.
    """
    input:
        table = segmentation_dir + '{well}_grid/{segmentation_type}.csv',
        composite = stitching_dir + '{well}_grid{grid_size}/grid_composite.json',
    output:
        table = segmentation_dir + '{well}_grid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/{segmentation_type}.csv'
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 2
    run:
        import pandas
        import numpy as np
        import constitch

        table = pandas.read_csv(input.table, index_col=0)
        composite = constitch.load(input.composite)
        grid_size, x, y = int(wildcards.grid_size), int(wildcards.x), int(wildcards.y)

        composite.boxes.positions[:,:2] *= phenotype_scale
        composite.boxes.positions[:,:2] //= bases_scale
        composite.boxes.sizes[:,:2] *= phenotype_scale
        composite.boxes.sizes[:,:2] //= bases_scale
        box = composite.boxes[x*grid_size+y]

        split_cells = 0

        contained = []
        for i, cell in table.iterrows():
            cellbox = constitch.BBox(point1=[cell.bbox_x1, cell.bbox_y1], point2=[cell.bbox_x2 // 2 * 2, cell.bbox_y2 // 2 * 2])
            # rounding down point2 of cellbox to avoid off by one error on the edge of the grid,
            # cause the grid will always be rounded down

            closest = np.linalg.norm(box.center - cellbox.center)
            closest_box = box
            closest_index = x * grid_size + y
            for j in range(len(composite.boxes)):
                dist = np.linalg.norm(composite.boxes[j].center - cellbox.center)
                if dist <= closest:
                    closest, closest_box = dist, composite.boxes[j]
                    closest_index = j

            #debug ('----- here ----', closest_index, closest_index//grid_size, closest_index%grid_size)

            if not closest_box.contains(cellbox):
                split_cells += 1
                debug ('split cell:', closest_index, closest_index//grid_size, closest_index%grid_size)
                #debug (closest_box.point1, closest_box.point2, cellbox.point1, cellbox.point2)
                #debug ('   ', np.linalg.norm(closest_box.center - cellbox.center))
                #for curbox in composite.boxes:
                    #debug (curbox.point1, curbox.point2)
                    #debug ('   ', cellbox.point1 - curbox.point1, curbox.point2 - cellbox.point2)
                    #debug ('   ', np.linalg.norm(cellbox.center - curbox.center))
            #assert closest_box.contains(cellbox)

            contained.append(closest_box is box)

            #is_contained = box.contains(cellbox)
            #for j in range(x*grid_size+y):
                #is_contained = is_contained and not composite.boxes[j].contains(cellbox)
            #contained.append(is_contained)

        assert split_cells < max(50, 0.2 * len(table.index)), (
                "{} out of {} cells in the well were split by the grid. "
                "Occasional cells are split if they are larger than overlap, but if "
                "this is too many consider increasing the overlap in config file".format(split_cells, len(table.index)))

        table = table[contained]
        table['bbox_x1'] -= box.position[0]
        table['bbox_x2'] -= box.position[0]
        table['bbox_y1'] -= box.position[1]
        table['bbox_y2'] -= box.position[1]
        table['xpos'] -= box.position[0]
        table['ypos'] -= box.position[1]
        table.to_csv(output.table)

def get_grid_filenames(wildcards):
    grid_size = int(wildcards.grid_size)
    numbers = ['{:02}'.format(i) for i in range(segmentation_grid_size)]
    return expand(segmentation_dir + '{well}_cellgrid' + str(segmentation_grid_size)
                + '/tile{x}x{y}y/{segmentation_type}_mask.tif', x=numbers, y=numbers, allow_missing=True)

def get_cells_mapping2(wildcards):
    if config['segmentation'].get('match_masks', False):
        return (segmentation_dir + '{well}_cellgrid' + str(segmentation_grid_size) + '/'
                + wildcards.segmentation_type.replace('cells', 'nuclei')
                + '_mappings.csv')
    return segmentation_dir + '{well}_cellgrid{grid_size,\d+}/{segmentation_type}_mappings.csv'

rule stitch_tile_segmentation:
    input:
        images = get_grid_filenames,
        composite = stitching_dir + '{well}_grid' + str(segmentation_grid_size) + '/grid_composite.json',
        #mappings = segmentation_dir + '{well}_cellgrid' + str(segmentation_grid_size) + '/{segmentation_type}_mappings.csv',
        mappings = get_cells_mapping2,
        composite2 = stitching_dir + '{well}_grid{grid_size}/grid_composite.json',
        table = segmentation_dir + '{well}_grid{grid_size}/tile{x}x{y}y/{segmentation_type}.csv',
    output:
        image = segmentation_dir + '{well}_grid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/{segmentation_type}_mask.tif',
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 2
    run:
        import tifffile
        import constitch

        composite = constitch.load(input.composite2)
        # scaling from base images to phenotype
        composite.boxes.positions[:,:2] *= phenotype_scale
        composite.boxes.positions[:,:2] //= bases_scale
        composite.boxes.sizes[:,:2] *= phenotype_scale
        composite.boxes.sizes[:,:2] //= bases_scale
        tile_index = int(wildcards.x) * int(wildcards.grid_size) + int(wildcards.y)

        tifffile.imwrite(output.image, stitch_segmentation_section(input.images,
                input.composite, input.mappings, composite.boxes[tile_index], input.table))

'''
rule split_grid_segmentation:
    """ Splits the segmentation masks into a smaller tile in a grid. Using the cell assignments made
    in the rule split_grid_table, only cells in the cell info table are copied into the segmentation mask
    for this tile. This ensures cells are not duplicated between tiles. As explained, overlap between
    tiles should be large enough that nearly all cells are fully contained in the tile they are assigned to.
    """
    input:
        image = segmentation_dir + '{well}_grid/{segmentation_type}.tif',
        composite = stitching_dir + '{well}_grid{grid_size}/grid_composite.json',
        table = segmentation_dir + '{well}_grid{grid_size}/tile{x}x{y}y/cells.csv',
    output:
        image = temp(segmentation_dir + '{well}_grid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/{segmentation_type}.tif'),
    wildcard_constraints:
        segmentation_type = '(cells|nuclei|cellsbases|nucleibases)_mask(_downscaled)?',
    run:
        import numpy as np
        import tifffile
        import pandas
        import constitch
        import skimage.measure

        image = tifffile.memmap(input.image, mode='r')
        table = pandas.read_csv(input.table, index_col=0)

        composite = constitch.load(input.composite)
        grid_size, x, y = int(wildcards.grid_size), int(wildcards.x), int(wildcards.y)

        downscaled = wildcards.segmentation_type.endswith('_downscaled')

        box = composite.boxes[x*grid_size+y]
        if not downscaled:
            box.position *= phenotype_scale
            box.position //= bases_scale
            box.size *= phenotype_scale
            box.size //= bases_scale

        debug (box)

        section = image[...,box.point1[0]:box.point2[0],box.point1[1]:box.point2[1]]
        #props = skimage.measure.regionprops(section)
        #props = {prop.label: prop for prop in props}

        dtype = [dtype for dtype in [np.uint16, np.uint32, np.uint64] if np.iinfo(dtype).max > len(table.index) + 1][0]
        #dtype = np.uint16 if np.iinfo(np.uint16).max > len(table.index) + 1 else np.uint32
        newimage = np.zeros(section.shape, dtype)

        if (wildcards.segmentation_type.count('cells') > 0
                or (wildcards.segmentation_type.count('nuclei') > 0
                and config['segmentation'].get('match_masks', False))):
            for newindex, (index, cell) in enumerate(table.iterrows()):
                x1, y1, x2, y2 = int(cell.bbox_x1), int(cell.bbox_y1), int(cell.bbox_x2), int(cell.bbox_y2)
                if downscaled:
                    x1, y1 = x1 * bases_scale // phenotype_scale, y1 * bases_scale // phenotype_scale
                    x2, y2 = x2 * bases_scale // phenotype_scale + 1, y2 * bases_scale // phenotype_scale + 1
                    # adding 1 to prevent rounding down cutting off segmentation
                x1, y1, x2, y2 = max(0, x1), max(0, y1), max(0, x2), max(0, y2)

                debug (x1, x2, y1, y2, box.point1, box.point2)
                mask = section[x1:x2,y1:y2] == index
                debug (x1, x2, y1, y2, mask.max(), mask.sum(), index, np.sum(section == index))
                #debug (' ', props[index].bbox)
                newimage[x1:x2,y1:y2][mask] = newindex + 1
        else:
            newimage, mapping, reverse_mapping = skimage.segmentation.relabel_sequential(section)
            newimage = newimage.astype(dtype)

        tifffile.imwrite(output.image, newimage)
'''

