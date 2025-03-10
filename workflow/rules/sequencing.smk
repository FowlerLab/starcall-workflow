import os
import glob
import re


rule segment_nuclei:
    input:
        sequencing_input_dir + '{prefix}/corrected_pt.tif',
        #sequencing_input_dir + '{prefix}/cycle' + phenotype_cycle + '.tif',
    output:
        sequencing_dir + '{prefix}/nuclei_mask.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 64 + 10000,
    threads: 2
    run:
        import numpy as np
        import tifffile
        import fisseq.segmentation

        data = tifffile.memmap(input[0])
        if data.shape[3] < 32:
            data = data.transpose(3,0,1,2)
        data = data.reshape(-1, *data.shape[2:])

        dapi = data[0]
        if np.all(dapi == 0):
            tifffile.imwrite(output[0], data[0])
        else:
            del data
            nuclei = fisseq.segmentation.segment_nuclei(dapi)
            debug ('Found', nuclei.max(), 'nuclei')
            tifffile.imwrite(output[0], nuclei)


rule segment_cells:
    input:
        sequencing_input_dir + '{prefix}/corrected_pt.tif',
        #sequencing_input_dir + '{prefix}/cycle' + phenotype_cycle + '.tif'
    output:
        #sequencing_output_dir + '{prefix}/nuclei_mask.tif',
        sequencing_dir + '{prefix}/cells_mask.tif',
        #sequencing_dir + '{prefix}/cells_{diam,\d+}diam_mask.tif',
        #sequencing_dir + '{prefix}/cells_{channel,\d+}channel_mask.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 16 + 10000,
        #cuda = 1,
    threads: 2
    run:
        import numpy as np
        import fisseq.segmentation
        import tifffile
        import logging
        import skimage.segmentation

        logging.basicConfig(level=logging.INFO)

        data = tifffile.memmap(input[0])
        debug (data.shape)
        if data.shape[3] < 32:
            data = data.transpose(3,0,1,2)
        debug (data.shape)
        data = data.reshape(-1, *data.shape[2:])

        debug(data.shape)

        dapi = data[0]
        cyto = data[2]
        if np.all(dapi == 0) or np.all(cyto == 0):
            tifffile.imwrite(output[0], data[0])
        else:
            #cyto = fisseq.segmentation.estimate_cyto(data[2:])
            del data
            #del full_well

            #cells = skimage.segmentation.expand_labels(nuclei, distance=cellpose_diameter / 3)
            #cells = fisseq.segmentation.segment_cyto_cellpose(cyto, dapi, diameter=cellpose_diameter * phenotype_scale)
            #diameter = int(wildcards.diam)
            cells = fisseq.segmentation.segment_cyto_cellpose(cyto, dapi, diameter=cellpose_diameter * phenotype_scale)
            #cells = fisseq.segmentation.segment_cyto_cellpose(cyto, dapi, diameter=cellpose_diameter * phenotype_scale, gpu=True)
            #nuclei = fisseq.segmentation.segment_nuclei(dapi)
            #cells = skimage.segmentation.expand_labels(cells, distance=5)

            debug ('Found', cells.max(), 'cells')
            #cells, nuclei = fisseq.segmentation.match_segmentations(cells, nuclei)
            #debug(f'found {cells.max()} nuclei/cells after reconciling')

            #tifffile.imwrite(output[0], nuclei)#, compression='deflate')
            tifffile.imwrite(output[0], cells)#, compression='deflate')

rule downscale_segmentation:
    input:
        sequencing_output_dir + '{prefix}/{mask}_mask.tif',
    output:
        sequencing_output_dir + '{prefix}/{mask,[^_]*}_mask_downscaled.tif',
    run:
        import tifffile
        import skimage.transform

        mask = tifffile.imread(input[0])

        tifffile.imwrite(output[0], skimage.transform.rescale(mask, 1/phenotype_scale, order=0))

ruleorder: downscale_segmentation > merge_grid
ruleorder: downscale_segmentation > segment_cells_bases

rule segment_cells_bases:
    input:
        #sequencing_input_dir + '{prefix}/corrected.tif'
        sequencing_input_dir + '{prefix}/raw.tif'
    output:
        sequencing_output_dir + '{prefix}/nuclei.tif',
        sequencing_output_dir + '{prefix}/cells.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 1.4 + 5000
    threads: 8
    run:
        import numpy as np
        import fisseq.segmentation
        import tifffile

        full_well = tifffile.memmap(input[0], mode='r')
        data = full_well[cycles.index(cellpose_cycle)].astype(np.float32)

        if np.all(data == 0):
            tifffile.imwrite(output[0], data[0])
            tifffile.imwrite(output[1], data[0])
        else:
            dapi = data[0]
            #cyto = data[2]
            cyto = fisseq.segmentation.estimate_cyto(data[2:])
            del data
            del full_well

            nuclei = fisseq.segmentation.segment_nuclei(dapi)
            cells = fisseq.segmentation.segment_cyto_cellpose(cyto, dapi, diameter=cellpose_diameter)
            #cells = skimage.segmentation.expand_labels(nuclei, cellpose_diameter // 2)

            debug(f'found {cells.max()} {nuclei.max()} nuclei/cells ')
            cells, nuclei = fisseq.segmentation.match_segmentations(cells, nuclei)
            debug(f'found {cells.max()} nuclei/cells after reconciling')

            tifffile.imwrite(output[0], nuclei)#, compression='deflate')
            tifffile.imwrite(output[1], cells)#, compression='deflate')


rule make_cell_overlay:
    input:
        image = sequencing_input_dir + '{prefix}/corrected_pt.tif',
        cells = sequencing_dir + '{prefix}/{segmentation_type}_mask.tif',
    output:
        qc_dir + '{prefix}/{segmentation_type}_overlay.tif',
        qc_dir + '{prefix}/{segmentation_type}_overlay.png',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 15 + 10000
    run:
        import numpy as np
        import tifffile
        import fisseq

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
        rgbimage = fisseq.utils.to_rgb8(image)
        tifffile.imwrite(output[1], rgbimage)


rule match_masks:
    input:
        cells = sequencing_dir + '{prefix}/cells_mask.tif',
        nuclei = sequencing_dir + '{prefix}/nuclei_mask.tif',
    output:
        cells = sequencing_output_dir + '{prefix}/cells_mask.tif',
        nuclei = sequencing_output_dir + '{prefix}/nuclei_mask.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 5 + 5000
    run:
        import collections
        import tifffile
        import numpy as np
        import skimage.measure
        import fisseq.segmentation

        cells = tifffile.imread(input.cells)
        nuclei = tifffile.imread(input.nuclei)

        cells, nuclei = fisseq.segmentation.match_segmentations(cells, nuclei)

        tifffile.imwrite(output.cells, cells)
        tifffile.imwrite(output.nuclei, nuclei)


rule make_cell_images:
    input:
        #image = sequencing_input_dir + '{prefix}/cycle' + phenotype_cycle + '.tif',
        image = sequencing_input_dir + '{prefix}/raw_pt.tif',
        cells = sequencing_output_dir + '{prefix}/cells_mask.tif',
        nuclei = sequencing_output_dir + '{prefix}/nuclei_mask.tif',
        #cells_nuclei = processing_dir + '{prefix}/cells_nuclei_mask.tif',
        cell_table = sequencing_output_dir + '{prefix}/cells.csv',
    output:
        cell_images = sequencing_output_dir + '{prefix}/cell_images_{window,\d+}.tif',
        mask_images = sequencing_output_dir + '{prefix}/mask_images_{window,\d+}.tif',
        #mask_nuclei_images = processing_dir + '{prefix}/mask_nuclei_images_{window,\d+}.tif',
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

"""
rule link_testing_prefix:
    input:
        sequencing_output_dir + '{prefix}/{file}',
    output:
        sequencing_output_dir + '{prefix}_test_{test,[^_]+}/{file}'
    wildcard_constraints:
        file = 'cells_mask.tif|cells_mask_downscaled.tif|cells.csv'
    shell:
        'cp -l {input[0]} {output[0]}'

ruleorder: tabulate_cells > link_testing_prefix
"""

rule tmp_link_file:
    input:
        sequencing_output_dir + 'well1_seqgrid5/tile02x02y_test_{test}/cells_reads.csv'
    output:
        'results/our_{test}.cells_reads.csv'
    localrule: True
    shell:
        'cp -l {input[0]} {output[0]}'

rule find_dots:
    input:
        sequencing_input_dir + '{prefix}/raw.tif'
    output:
        sequencing_dir + '{prefix}/bases.csv'
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 10 + 15000
    threads: 4
    run:
        import numpy as np
        import tifffile
        import fisseq.dotdetection
        import fisseq.correction
        import skimage.morphology

        min_sigma = 1#2
        max_sigma = 2#3
        num_sigma = 7

        full_well = tifffile.memmap(input[0], mode='r')
        image = np.array(full_well[...,2:,:,:], dtype=np.float32)
        del full_well

        if np.all(image == 0):
            tifffile.imwrite(output.dots, image[0,0])
            tifffile.imwrite(output.dotvalues, image[0,0])
        else:

            print_mem('extract_bases', resources.mem_mb)
            debug('image', image.shape)

            #dot_filter = fisseq.dotdetection.dot_filter(image, copy=False)
            #tifffile.imwrite('tmp_dots_dot_filter.tif', dot_filter)

            """
            oldshape = dot_filter.shape
            dot_filter = dot_filter.reshape(oldshape[0] * oldshape[1], oldshape[2], oldshape[3])
            newvalues = np.empty_like(dot_filter)
            mask = skimage.morphology.disk(1)
            for i in range(dot_filter.shape[0]):
                skimage.morphology.dilation(dot_filter[i], mask, out=newvalues[i])
            dot_filter = newvalues.reshape(oldshape)
            #"""

            poses, values = fisseq.dotdetection.detect_dots3(
                image,
                min_sigma = min_sigma,
                max_sigma = max_sigma,
                num_sigma = num_sigma,
                #copy = False,
            )

            #debug('dot_filter', dot_filter.shape)

            # take read values from original image or filtered image
            #values = image[:,:,poses[:,0],poses[:,1]]
            #values = dot_filter[:,:,poses[:,0],poses[:,1]]

            values = values.transpose(1,0,2)
            debug("correcting", values.shape)
            dye_matrix = fisseq.correction.estimate_dye_matrix(values)
            corrected = fisseq.correction.color_correct(values)
            debug("corrected")
            fisseq.correction.crosstalk_plot(values, corrected, dye_matrix, name=wildcards.prefix.replace('/','_'))
            #values = corrected
            values = values.transpose(1,0,2)

            values = values.reshape(image.shape[0] * image.shape[1], -1).T
            np.savetxt(output[0], np.concatenate((poses, values), axis=1), delimiter=',', fmt='%f')



ruleorder: merge_grid > find_dots
ruleorder: merge_grid > segment_cells
#ruleorder: link_grid > find_dots
#ruleorder: link_grid > segment_cells
ruleorder: merge_grid > segment_cells_bases
ruleorder: segment_cells > segment_cells_bases


rule annotate_dots:
    input:
        #input_dir + '{prefix}/cycle' + cycles[0] + '.tif',
        #sequencing_input_dir + '{prefix}/corrected.tif',
        sequencing_input_dir + '{prefix}/raw.tif',
        sequencing_dir + '{prefix}/bases.csv',
    output:
        qc_dir + '{prefix}/annotated.tif',
    resources:
        mem_mb = 25000
    run:
        import numpy as np
        import fisseq.utils
        import tifffile

        #image = tifffile.imread(input[0])
        image = tifffile.memmap(input[0], mode='r')[0]
        poses = np.loadtxt(input[1], delimiter=',')[:,:2].astype(int)
        tifffile.imwrite(output[0], fisseq.utils.mark_dots(image, poses))


rule tabulate_cells:
    input:
        cells = sequencing_output_dir + '{prefix}/{segmentation_type}_mask.tif',
    output:
        table = sequencing_output_dir + '{prefix}/{segmentation_type}.csv',
    wildcard_constraints:
        segmentation_type = 'cells|nuclei'
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 1.5
    run:
        import numpy as np
        import tifffile
        import skimage.measure
        import fisseq.sequencing
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

ruleorder: tabulate_cells > merge_grid

rule call_reads:
    input:
        cell_table = sequencing_output_dir + '{prefix}/{segmentation_type}.csv',
        #bases = sequencing_dir + '{prefix}/bases.csv',
        bases = sequencing_dir + '{prefix}/bases.csv',
        cells = sequencing_output_dir + '{prefix}/{segmentation_type}_mask_downscaled.tif',
    output:
        sequencing_dir + '{prefix}/{segmentation_type}_reads_partial.csv',
        #sequencing_dir + '{prefix}/{segmentation_type}_reads_partial.csv',
    wildcard_constraints:
        segmentation_type = 'cells|nuclei'
    resources:
        mem_mb = lambda wildcards, input: 10000 + input.size_mb * 5
    run:
        import numpy as np
        import tifffile
        import skimage.measure
        import fisseq.sequencing
        import pandas

        dist_thresh = 1.5#float(wildcards.dist_limit)

        poses = np.loadtxt(input.bases, delimiter=',')
        poses, values = poses[:,:2].astype(int), poses[:,2:]
        values = values.reshape(values.shape[0], -1, 4)

        num_bases = values.shape[1]
        cells = tifffile.imread(input.cells)
        cell_table = pandas.read_csv(input.cell_table, index_col=0)

        reads, counts = fisseq.sequencing.call_reads(cells, poses, values, distance_threshold=dist_thresh)

        table = {}
        index = cell_table.index
        reads = [reads[i] for i in index]
        counts = [counts[i] for i in index]

        fisseq.utils.write_multicolumn(table, 'read', reads, length_column='read_count')
        fisseq.utils.write_multicolumn(table, 'count', counts, length_column='read_count')

        table['total_count'] = list(map(sum, counts))

        table = pandas.DataFrame(table, index=index)
        table.to_csv(output[0])



def get_aux_data(wildcards, prefix=None):
    prefix = wildcards.prefix if prefix is None else prefix

    if prefix != '': prefix = prefix + '.'

    files = []
    for base_dir in (sequencing_dir, input_dir):
        pattern = base_dir + '{prefix}/{segmentation_type}.auxdata/*.csv'.format(prefix=prefix, segmentation_type=wildcards.segmentation_type)
        files.extend(sorted(glob.glob(pattern)))
        pattern = base_dir + '{prefix}/auxdata/*.csv'.format(prefix=prefix, segmentation_type=wildcards.segmentation_type)
        files.extend(sorted(glob.glob(pattern)))

    #if re.fullmatch('(tile.+)|(well.+)|(cycle.+)', os.path.basename(prefix)):
    if prefix:
        files.extend(get_aux_data(wildcards, prefix=os.path.dirname(prefix)))

    return files

rule merge_final_tables:
    input:
        cell_table = sequencing_dir + '{prefix}/{segmentation_type}_reads_partial.csv',
        aux_data = get_aux_data,
    output:
        full_table = sequencing_output_dir + '{prefix}/{segmentation_type}_reads.csv',
    wildcard_constraints:
        segmentation_type = 'cells|nuclei'
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 250
    run:
        import pandas
        import numpy as np
        import fisseq.sequencing
        import warnings

        debug ('begin ')
        cell_table = pandas.read_csv(input.cell_table, index_col=0)
        debug (cell_table)
        cell_cols = cell_table.columns.copy()

        def join_barcode(cell_table, aux_table):
            barcodes = []
            for barc in aux_table.index:
                new_barcodes = barc.split('-')
                for new_barc in new_barcodes:
                    if len(new_barc) != len(cycles):
                        warnings.warn('Barcode not the right length {} {}'.format(len(new_barc), len(cycles)))
                barcodes.append([subbarc[:len(cycles)] for subbarc in new_barcodes])

            barcodes = np.array(barcodes)
            lengths = np.array(list(map(len, barcodes.flat)))
            debug (lengths)
            debug (lengths.min(), lengths.mean(), lengths.max())

            library = fisseq.sequencing.BarcodeLibrary(barcodes)

            reads = []
            counts = []
            index = 0

            while 'read_{}'.format(index) in cell_table.columns and cell_table['count_{}'.format(index)].sum() > 0:
                reads.append(list(cell_table['read_{}'.format(index)].to_numpy()))
                counts.append(cell_table['count_{}'.format(index)].to_numpy())
                index += 1

            while 'barcode_{}'.format(index) in cell_table.columns and cell_table['count_{}'.format(index)].sum() > 0:
                reads.append(list(cell_table['barcode_{}'.format(index)].to_numpy()))
                counts.append(cell_table['count_{}'.format(index)].to_numpy())
                index += 1
            
            debug('combingin reads and stuff', len(reads), len(counts))
            reads = np.stack(reads, axis=-1)
            counts = np.stack(counts, axis=-1)
            lengths = np.array(list(map(len, reads.flat)))
            debug (reads.dtype)
            debug (barcodes.dtype)

            debug ('starting matching')
            #indices, edit_distances = fisseq.sequencing.match_barcodes(reads, counts, barcodes, max_edit_distance=99999, return_distances=True, debug=True, progress=True)
            edit_distances, indices = library.nearest(reads, counts)
            debug ('  done')
            debug (np.sum(indices != -1) / len(indices))

            new_rows = [pandas.Series(dtype=object) if i == -1 else aux_table.iloc[i,:] for i in indices]
            new_table = pandas.DataFrame(new_rows, index=cell_table.index)
            new_table['editDistance'] = edit_distances
            return cell_table.join(new_table)

        for path in input.aux_data:
            aux_table = pandas.read_csv(path)
            debug (aux_table)

            if path.endswith(wildcards.segmentation_type + '.csv'):
                cell_table = cell_table.join(aux_table.set_index(aux_table.columns[0]), how='left')
            elif path.endswith('reads.csv') or path.endswith('barcodes.csv'):
                debug('staring join')
                cell_table = join_barcode(cell_table, aux_table.set_index(aux_table.columns[0]))
                #cell_table = cell_table.reset_index().merge(aux_table, how='left', left_on='barcode_0', right_on=aux_table.columns[0]).set_index('index')
                #cell_table = cell_table.join(aux_table.set_index(aux_table.columns[0]), how='left', on='barcode_0')

            elif len(aux_table.index) == 1:
                for col in aux_table.columns:
                    cell_table[col] = [aux_table[col][0]] * len(cell_table.index)
            else:
                shared_cols = list(set(cell_cols) & set(aux_table.columns))
                debug (shared_cols)
                if len(shared_cols) != 0:
                    cell_table = cell_table.join(aux_table.set_index(shared_cols), how='left', on=shared_cols)
                else:

                    cell_cols = [col for col in aux_table.columns if col.lower() in ['cell', 'cellid', 'cellindex', 'cell_id', 'cell_index']]
                    barc_cols = [col for col in aux_table.columns if col.lower() in ['barcode', 'read']]

                    if len(cell_cols) != 0:
                        cell_table = cell_table.join(aux_table.set_index(cell_cols[0]), how='left')
                    elif len(barc_cols) != 0:
                        cell_table = join_barcode(cell_table, aux_table.set_index(barc_cols))

            debug (cell_table)

        cell_table.to_csv(output.full_table)

ruleorder: call_reads > merge_grid
ruleorder: match_masks > merge_grid
ruleorder: merge_final_tables > merge_grid
ruleorder: calc_features > merge_grid

