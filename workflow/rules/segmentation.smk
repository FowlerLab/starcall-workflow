import os
import glob
import re


def get_segmentation_pt(wildcards):
    path = wildcards.path_nogrid.replace('_cellgrid', '_grid')
    if config['segmentation']['use_corrected']:
        return stitching_dir + path + '/corrected_pt.tif'
    else:
        return stitching_dir + path + '/raw_pt.tif'

rule segment_nuclei:
    input:
        get_segmentation_pt,
    output:
        segmentation_dir + '{path_nogrid}/nuclei{nuclearchannel}_mask_unmatched.tif',
    params:
        nuclearchannel = parse_param('nuclearchannel', config['segmentation']['channels'][0])
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

        nuclearchannel = channel_index(params.nuclearchannel, kind='phenotyping')

        data = tifffile.memmap(input[0], mode='r')
        if data.shape[3] < 32:
            data = data.transpose(3,0,1,2)
        data = data.reshape(-1, *data.shape[2:])

        dapi = data[nuclearchannel]
        if np.all(dapi == 0):
            tifffile.imwrite(output[0], data[0])
        else:
            del data
            nuclei = starcall.segmentation.segment_nuclei(dapi)
            debug ('Found', nuclei.max(), 'nuclei')
            tifffile.imwrite(output[0], nuclei)


rule segment_cells:
    input:
        get_segmentation_pt,
    output:
        segmentation_dir + '{path_nogrid}/cells{diameter}{nuclearchannel}{cytochannel}_mask_unmatched.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 16 + 6000,
        cuda = 1,
    params:
        diameter = parse_param('diameter', config['segmentation']['diameter']),
        nuclearchannel = parse_param('nuclearchannel', config['segmentation']['channels'][0]),
        cytochannel = parse_param('cytochannel', config['segmentation']['channels'][1]),
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

        diameter = int(wildcards.diameter)
        nuclearchannel = channel_index(params.nuclearchannel, kind='phenotyping')
        cytochannel = channel_index(params.cytochannel, kind='phenotyping')

        if type(nuclearchannel) == str:
            nuclearchannel = config['phenotyping_channels'].index(nuclearchannel)

        if type(cytochannel) == str:
            cytochannel = config['phenotyping_channels'].index(cytochannel)

        logging.basicConfig(level=logging.INFO)

        data = tifffile.memmap(input[0], mode='r')
        debug (data.shape)
        if data.shape[3] < 32:
            data = data.transpose(3,0,1,2)
        debug (data.shape)
        data = data.reshape(-1, *data.shape[2:])

        debug(data.shape)

        dapi = data[nuclearchannel]
        cyto = data[cytochannel]
        if np.all(dapi == 0) or np.all(cyto == 0):
            tifffile.imwrite(output[0], data[0])
        else:
            del data
            #del full_well

            cells = starcall.segmentation.segment_cyto_cellpose(
                cyto, dapi,
                diameter = config['segmentation']['diameter'],
                gpu=resources.cuda==1,
            )

            debug ('Found', cells.max(), 'cells')

            tifffile.imwrite(output[0], cells)#, compression='deflate')


def get_segmentation_bases(wildcards):
    path = wildcards.path_nogrid.replace('_cellgrid', '_grid')
    return stitching_dir + path + '/raw.tif'

rule segment_cells_bases:
    input:
        get_segmentation_bases,
    output:
        segmentation_dir + '{path_nogrid}/cellsbases{diameter}{nuclearchannel}_mask_unmatched.tif',
    params:
        diameter = parse_param('diameter', config['segmentation']['diameter']),
        nuclearchannel = parse_param('nuclearchannel', config['segmentation']['channels'][0])
    wildcard_constraints:
        diameter = '|_diameter\d+',
        nuclearchannel = '|_nuclearchannel' + sequencing_channel_regex,
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 8 + 10000,
        cuda = 1,
    threads: 8
    run:
        import numpy as np
        import starcall.segmentation
        import tifffile

        diameter = params.diameter
        nuclearchannel = channel_index(params.nuclearchannel, kind='sequencing')

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

            cells = starcall.segmentation.segment_cyto_cellpose(
                cyto, dapi,
                diameter = config['segmentation']['diameter'] * bases_scale // phenotype_scale,
                gpu=resources.cuda==1,
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
        mem_mb = lambda wildcards, input: input.size_mb * 8 + 10000
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
    input:
        segmentation_dir + '{path_nogrid}/{segmentation_type}_mask.tif',
    output:
        segmentation_dir + '{path_nogrid}/{segmentation_type}_mask_downscaled.tif',
    run:
        import tifffile
        import skimage.transform

        mask = tifffile.imread(input[0])

        if wildcards.segmentation_type in ('cellsbases', 'nucleibases'):
            tifffile.imwrite(output[0], mask)
        else:
            tifffile.imwrite(output[0], skimage.transform.rescale(mask, bases_scale/phenotype_scale, order=0))


rule make_cell_overlay:
    input:
        image = get_segmentation_pt,
        cells = segmentation_dir + '{path_nogrid}/{segmentation_type}_mask{params}.tif',
    output:
        qc_dir + '{path_nogrid}/{segmentation_type}_overlay{params}.tif',
        qc_dir + '{path_nogrid}/{segmentation_type}_overlay{params}.png',
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


rule make_cell_images:
    input:
        #image = stitching_dir + '{path}/cycle' + phenotype_cycle + '.tif',
        image = get_segmentation_pt,
        cells = segmentation_dir + '{path_nogrid}/cells_mask.tif',
        nuclei = segmentation_dir + '{path_nogrid}/nuclei_mask.tif',
        #cells_nuclei = processing_dir + '{path}/cells_nuclei_mask.tif',
        cell_table = segmentation_dir + '{path_nogrid}/cells.csv',
    output:
        cell_images = segmentation_dir + '{path_nogrid}/cell_images_{window,\d+}.tif',
        mask_images = segmentation_dir + '{path_nogrid}/mask_images_{window,\d+}.tif',
        #mask_nuclei_images = processing_dir + '{path}/mask_nuclei_images_{window,\d+}.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 2.5 + 5000
    run:
        import numpy as np
        import tifffile
        import pandas
        
        #cell_table = np.genfromtxt(input.cell_table, delimiter=',', dtype=None, names=None)
        cell_table = pandas.read_csv(input.cell_table, index_col=0)
        cells = tifffile.imread(input.cells)
        nuclei = tifffile.imread(input.nuclei)
        image = tifffile.imread(input.image)

        image = image.reshape(-1, *image.shape[2:])
        debug (image.shape)

        window = int(wildcards.window)
        window_low = window // 2
        window_high = window - window_low

        cell_images = np.zeros((len(cell_table), image.shape[0], window, window), image.dtype)
        mask_images = np.zeros((len(cell_table), 2, window, window), dtype=bool)

        for i, cell_index in enumerate(cell_table.index):
            debug (cell_index)
            centroid = int(cell_table['xpos'][cell_index]) * 2, int(cell_table['ypos'][cell_index]) * 2
            x1, x2, y1, y2 = centroid[0] - window_low, centroid[0] + window_high, centroid[1] - window_low, centroid[1] + window_high
            x1, x2, y1, y2 = max(0, x1), min(image.shape[1], x2), max(0, y1), min(image.shape[2], y2)
            subset = image[:,x1:x2,y1:y2]
            mask = cells[x1:x2,y1:y2] == cell_index, nuclei[x1:x2,y1:y2] == cell_index
            x1, x2 = window_low - (centroid[0] - x1), window_low + (x2 - centroid[0])
            y1, y2 = window_low - (centroid[1] - y1), window_low + (y2 - centroid[1])
            cell_images[i,:,x1:x2,y1:y2] = subset
            mask_images[i,0,x1:x2,y1:y2] = mask[0]
            mask_images[i,1,x1:x2,y1:y2] = mask[1]

        tifffile.imwrite(output.cell_images, cell_images)
        tifffile.imwrite(output.mask_images, mask_images)


rule tabulate_cells:
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

rule merge_grid_segmentation:
    input:
        images = get_grid_filenames,
        composite = stitching_dir + '{well}_grid{grid_size}/grid_composite.json',
    output:
        image = segmentation_dir + '{well}_cellgrid{grid_size,\d+}/{segmentation_type}_mask_unmatched.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 10 + 25000
    run:
        import numpy as np
        import tifffile
        import constitch
        import pandas

        composite = constitch.load(input.composite)

        for i,path in enumerate(input.images):
            composite.images[i] = tifffile.imread(path)

        # scaling from base images to phenotype
        composite.boxes.positions[:,:2] *= phenotype_scale
        composite.boxes.positions[:,:2] //= bases_scale
        composite.boxes.sizes[:,:2] *= phenotype_scale
        composite.boxes.sizes[:,:2] //= bases_scale
        merger = constitch.MaskMerger()

        full_image = composite.stitch(merger=merger)
        del composite
        tifffile.imwrite(output.image, full_image)


segmentation_grid_size = config.get('segmentation_grid_size', 1)

def get_grid_size_file(wildcards):
    grid_size = config.get('segmentation_{}_grid_size'.format(wildcards.segmentation_type), segmentation_grid_size)
    if grid_size == 1:
        return segmentation_dir + '{well}/{segmentation_type}_mask_unmatched.tif'
    return segmentation_dir + '{well}_cellgrid' + str(grid_size) + '/{segmentation_type}_mask_unmatched.tif',

rule link_merged_grid:
    input:
        get_grid_size_file,
    output:
        segmentation_dir + '{well}_grid/{segmentation_type}_mask_unmatched.tif',
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"

ruleorder: link_merged_grid > segment_cells
ruleorder: link_merged_grid > segment_nuclei


rule split_grid_table:
    input:
        table = segmentation_dir + '{well}_grid/cells.csv',
        composite = stitching_dir + '{well}_grid{grid_size}/grid_composite.json',
    output:
        table = segmentation_dir + '{well}_grid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/cells.csv'
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

        contained = []
        for i, cell in table.iterrows():
            cellbox = constitch.BBox(point1=[cell.bbox_x1, cell.bbox_y1], point2=[cell.bbox_x2, cell.bbox_y2])

            closest = np.linalg.norm(box.center - cellbox.center)
            closest_box = box
            #if i == 5487:
                #debug ('begin', box, cellbox, box.contains(cellbox), closest)
            for j in range(len(composite.boxes)):
                dist = np.linalg.norm(composite.boxes[j].center - cellbox.center)
                #if i == 5487:
                    #debug (composite.boxes[j], cellbox, composite.boxes[j].contains(cellbox), dist)
                if dist <= closest:
                    closest, closest_box = dist, composite.boxes[j]

            assert closest_box.contains(cellbox)

            contained.append(closest_box is box)

            #is_contained = box.contains(cellbox)
            #for j in range(x*grid_size+y):
                #is_contained = is_contained and not composite.boxes[j].contains(cellbox)
            #contained.append(is_contained)

        table = table[contained]
        table['bbox_x1'] -= box.position[0]
        table['bbox_x2'] -= box.position[0]
        table['bbox_y1'] -= box.position[1]
        table['bbox_y2'] -= box.position[1]
        table['xpos'] -= box.position[0]
        table['ypos'] -= box.position[1]
        table.to_csv(output.table)


rule split_grid_segmentation:
    input:
        image = segmentation_dir + '{well}_grid/{segmentation_type}.tif',
        composite = stitching_dir + '{well}_grid{grid_size}/grid_composite.json',
        table = segmentation_dir + '{well}_grid{grid_size}/tile{x}x{y}y/cells.csv',
    output:
        image = segmentation_dir + '{well}_grid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/{segmentation_type}.tif',
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

        dtype = [dtype for dtype in [np.uint8, np.uint16, np.uint32] if np.iinfo(dtype).max > len(table.index) + 1][0]
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

                debug (x1, x2, y1, y2, box.point1, box.point2)
                mask = section[x1:x2,y1:y2] == index
                debug (x1, x2, y1, y2, mask.max(), mask.sum(), index, np.sum(section == index))
                #debug (' ', props[index].bbox)
                newimage[x1:x2,y1:y2][mask] = newindex + 1
        else:
            newimage, mapping, reverse_mapping = skimage.segmentation.relabel_sequential(section)
            newimage = newimage.astype(dtype)

        tifffile.imwrite(output.image, newimage)


