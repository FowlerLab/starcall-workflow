import os
import sys
import glob

##################################################
## Phenotyping
##################################################

def get_phenotyping_pt(wildcards):
    if config['phenotyping']['use_corrected']:
        return phenotyping_dir + '{path}/corrected_pt.tif'
    else:
        return phenotyping_dir + '{path}/raw_pt.tif'


rule make_cell_images:
    """ Create small crops of each cell segmented, useful if the cell images are being provided
    to an image embedding network or similar procedure. The images are cropped to a set size, specified
    in the output filename, with the centroid of the cell in the center of the image. The output
    tifffile has shape (num_cells, num_channels, window_size, window_size). In addition the cell masks
    are cropped and saved as boolean masks with shape (num_cells, window_size, window_size)
    """
    input:
        image = get_phenotyping_pt,
        cells = phenotyping_dir + '{path}/{segmentation_type}_mask.tif',
        cell_table = phenotyping_dir + '{path}/{segmentation_type}.csv',
    output:
        cell_images = phenotyping_dir + '{path}/{segmentation_type}_crops_{window,\d+}.tif',
        mask_images = phenotyping_dir + '{path}/{segmentation_type}_mask_crops_{window,\d+}.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 2.5 + 5000
    run:
        import numpy as np
        import tifffile
        import pandas
        
        #cell_table = np.genfromtxt(input.cell_table, delimiter=',', dtype=None, names=None)
        cell_table = pandas.read_csv(input.cell_table, index_col=0)

        if len(cell_table.index) == 0:
            os.system('touch {}'.format(output.cell_images))
            os.system('touch {}'.format(output.mask_images))

        else:
            cells = tifffile.imread(input.cells)
            image = tifffile.imread(input.image)

            image = image.reshape(-1, *image.shape[2:])
            debug (image.shape)

            window = int(wildcards.window)
            window_low = window // 2
            window_high = window - window_low

            cell_images = np.zeros((len(cell_table), image.shape[0], window, window), image.dtype)
            mask_images = np.zeros((len(cell_table), window, window), dtype=np.uint8) # bool images are not memory mappable

            for i, cell_index in enumerate(cell_table.index):
                debug (cell_index)
                centroid = int(cell_table['xpos'][cell_index]), int(cell_table['ypos'][cell_index])
                x1, x2, y1, y2 = centroid[0] - window_low, centroid[0] + window_high, centroid[1] - window_low, centroid[1] + window_high
                x1, x2, y1, y2 = max(0, x1), min(cells.shape[0], x2), max(0, y1), min(cells.shape[1], y2)
                debug (x1, x2, y1, y2)
                subset = image[:,x1:x2,y1:y2]
                mask = cells[x1:x2,y1:y2] == i + 1
                x1, x2 = window_low - (centroid[0] - x1), window_low + (x2 - centroid[0])
                y1, y2 = window_low - (centroid[1] - y1), window_low + (y2 - centroid[1])
                cell_images[i,:,x1:x2,y1:y2] = subset
                mask_images[i,x1:x2,y1:y2] = mask

            tifffile.imwrite(output.cell_images, cell_images)
            tifffile.imwrite(output.mask_images, mask_images)

rule extract_embeddings:
    input:
        cell_images = phenotyping_dir + '{path}/cells_crops_100.tif',
        cells = phenotyping_dir + '{path}/cells.csv',
    output:
        embeddings = phenotyping_dir + '{path}/embeddings_morphem{cycle,|_cycle\d+}.csv',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 7 + 5000,
        cuda = 1,
    run:
        import starcall.embedding
        import pandas
        import tifffile

        images = tifffile.imread(input.cell_images)
        debug (images.shape)
        features = starcall.embedding.morphem(images, device='cuda', debug=True, progress=True)

        cells_table = pandas.read_csv(input.cells, index_col=0)

        feature_table = pandas.DataFrame({'morphem-{:04}'.format(i): features[:,i] for i in range(features.shape[1])}, index=cells_table.index)
        feature_table.to_csv(output.embeddings)

rule calc_features:
    input:
        cell_table = phenotyping_dir + '{path}/cells.csv',
        cells = phenotyping_dir + '{path}/cells_mask.tif',
        nuclei = phenotyping_dir + '{path}/nuclei_mask.tif',
        image = get_phenotyping_pt,
    output:
        features = phenotyping_dir + '{path}/features.csv'
    resources:
        mem_mb = lambda wildcards, input, attempt: input.size_mb * 2 + 5000
    run:
        import tifffile
        import numpy as np
        import pandas
        import skimage.measure
        import starcall.utils

        cell_table = pandas.read_csv(input.cell_table, index_col=0)
        cells, nuclei = tifffile.imread(input.cells), tifffile.imread(input.nuclei)
        image = tifffile.imread(input.image)
        image = image.reshape(-1, *image.shape[2:])

        all_props = skimage.measure.regionprops(cells, image.transpose(1,2,0))
        all_props = {props.label: props for props in all_props}

        features = {}

        for i, (cell_index, cell) in enumerate(starcall.utils.simple_progress(list(cell_table.iterrows()))):
            x1, y1 = int(cell.bbox_x1), int(cell.bbox_y1) + 1
            x2, y2 = int(cell.bbox_x2), int(cell.bbox_y2) + 1
            cell_mask = cells[x1:x2,y1:y2] == i + 1
            nucleus_mask = nuclei[x1:x2,y1:y2] == i + 1
            image_section = image[:,x1:x2,y1:y2]
            props = all_props[i+1]

            for prop in props:
                if type(props[prop]) in (float, int):
                    features.setdefault(prop, []).append(props[prop])

            for maskname, mask in [('cell', cell_mask), ('nucleus', nucleus_mask), ('cytoplasm', cell_mask & ~nucleus_mask)]:
                for channel in range(len(config['phenotyping_channels'])):
                    basename = '{}_ch{}'.format(maskname, channel)
                    masked_section = image_section[channel,mask]

                    percentiles = [0,1,5,50,95,99,100]
                    values = np.array([0] * len(percentiles))
                    area_mean, area_sum = 0, 0
                    if masked_section.size > 0:
                        values = np.percentile(masked_section, percentiles)
                        area_mean = masked_section.mean()
                        area_sum = masked_section.sum()

                    for percent, val in zip(percentiles, values):
                        feature_name = basename + '_{}percentile'.format(percent)
                        features.setdefault(feature_name, []).append(val)

                    features.setdefault(basename + '_min', []).append(values[0])
                    features.setdefault(basename + '_max', []).append(values[-1])
                    features.setdefault(basename + '_mean', []).append(area_mean)
                    features.setdefault(basename + '_sum', []).append(area_sum)

        features = pandas.DataFrame(features, index=cell_table.index)
        features.to_csv(output.features)

##################################################
## Phenotyping with cellprofiler
##################################################

rule extract_cellprofiler_channel:
    """ Copies a single channel from the phenotype images, as cellprofiler requires each channel to be as separate image
    """
    input:
        image = get_phenotyping_pt,
    output:
        image = temp(phenotyping_dir + '{path}/cellprofiler{cycle}/channel{channel,\d+\.\d+}.tif'),
    wildcard_constraints:
        cycle = '|cycle\d+',
    run:
        import numpy as np
        import tifffile

        image = tifffile.memmap(input.image, mode='r')
        cycle, chan = wildcards.channel.split('.')
        tifffile.imwrite(output.image, image[int(cycle),int(chan)])

max_num_channels = max(len(channels) for channels in config['phenotyping_channels'])

def get_channels(wildcards):
    params_channels = config['phenotyping'].get('channels', None)
    if wildcards.cycle != '':
        cycle = int(wildcards.cycle[5:])
        channel_indices = [(cycle, i) for i in range(len(config['phenotyping_channels'][cycle]))]
    elif params_channels is not None:
        channel_indices = [channel_index_phenotyping(chan) for chan in params_channels]
    else:
        channel_indices = []
        for cycle, channels in enumerate(config['phenotyping_channels'][0]):
            channel_indices.extend((cycle, i) for i in range(len(channels)))
    return [phenotyping_dir + '{path}/cellprofiler{cycle}/channel' + str(cycle) + '.' + str(chan) + '.tif' for cycle,chan in channel_indices]


rule copy_cellprofiler_files:
    """ Ensures all input for cellprofiler is in the correct format and location. writes the file list
    that is read in by the cellprofiler pipeline.
    """
    input:
        #image = stitching_dir + '{path}/cycle' + phenotype_cycle + '.tif',
        #image = get_phenotyping_pt,
        images = get_channels,
        cells = phenotyping_dir + '{path}/cells_mask.tif',
        nuclei = phenotyping_dir + '{path}/nuclei_mask.tif',
        #puncta = phenotyping_dir + '{path}/puncta_mask.tif',
        #lines = phenotyping_dir + '{path}/line_mask.tif',
    output:
        file_list = phenotyping_dir + '{path}/cellprofiler{cycle}/files.csv',
        #images = expand(phenotyping_dir + '{path}/cellprofiler{cycle}/channel{channel}.tif', channel=range(max_num_channels), allow_missing=True),
        cells = temp(phenotyping_dir + '{path}/cellprofiler{cycle}/cells.tif'),
        nuclei = temp(phenotyping_dir + '{path}/cellprofiler{cycle}/nuclei.tif'),
        #puncta = phenotyping_dir + '{path}/cellprofiler/puncta.tif',
        #lines = phenotyping_dir + '{path}/cellprofiler/lines.tif',
    wildcard_constraints:
        cycle = '|cycle\d+',
    params:
        #cycle = parse_param('cycle', None)
        channels = config['phenotyping'].get('channels', None),
    run:
        import numpy as np
        import tifffile

        with open(output.file_list, 'w') as ofile:
            if len(input) - len(input.images) > 2:
                ofile.write(','.join(['FileName_CH{}'.format(i) for i in range(len(input.images))]) + ',FileName_Cells,FileName_Nuclei,FileName_Puncta,FileName_Line\n')
            else:
                ofile.write(','.join(['FileName_CH{}'.format(i) for i in range(len(input.images))]) + ',FileName_Cells,FileName_Nuclei\n')

            for i, path in enumerate(input.images):
                ofile.write(os.path.basename(path) + ',')

            startindex = len(input.images)
            for i, path, outpath in zip(range(len(input) - startindex), input[startindex:], output[1:]):
                debug (i, path, outpath)
                with tifffile.TiffFile(path) as cells_file:
                    dtype = cells_file.pages[0].dtype

                if dtype != np.uint16:
                    image = tifffile.imread(path)
                    if image.max() > np.iinfo(np.uint16).max:
                        assert dtype == np.uint16, (
                            "Segmentation input to cellprofiler must be in uint16 form. "
                            "'{}' has dtype {}, this may be due to there being more than 65535 "
                            "cells in the image. To solve this increase the grid size for phenotyping "
                            " in config.yaml.".format(path, dtype))
                    else:
                        tifffile.imwrite(outpath, image.astype(np.uint16))
                else:
                    #os.symlink(os.path.relpath(path, os.path.dirname(outpath)), outpath)
                    os.link(path, outpath)

                if i != 0:
                    ofile.write(',')
                ofile.write(os.path.basename(outpath))

            ofile.write('\n')


def find_pipeline(wildcards):
    pipeline = glob.glob('*{}.cppipe'.format(wildcards.pipeline))
    if len(pipeline) != 1:
        return wildcards.pipeline + '.cppipe'
    return pipeline[0]

rule run_cellprofiler:
    input:
        file_list = phenotyping_dir + '{path}/cellprofiler{cycle}/files.csv',
        #images = expand(phenotyping_dir + '{path}/cellprofiler{cycle}/channel{channel}.tif', channel=range(len(config['phenotyping_channels'])), allow_missing=True),
        images = get_channels,
        cells = phenotyping_dir + '{path}/cellprofiler{cycle}/cells.tif',
        nuclei = phenotyping_dir + '{path}/cellprofiler{cycle}/nuclei.tif',
        #puncta = phenotyping_dir + '{path}/cellprofiler/puncta.tif',
        #lines = phenotyping_dir + '{path}/cellprofiler/lines.tif',
        pipeline = find_pipeline,
        #pipeline = '{pipeline}.cppipe',
    output:
        #data = phenotyping_dir + '{path}/cellprofiler_{pipeline,[^./]+}.csv',
        #mark = phenotyping_dir + '{path}/cellprofiler/{pipeline,[^./]+}/mark',
        cell_file = phenotyping_dir + '{path}/cellprofiler{cycle,|cycle\d+}/{pipeline}/Cells.csv'
    params:
        cellprofiler_executable = config['phenotyping'].get('cellprofiler_executable', 'cellprofiler'),
    resources:
        mem_mb = lambda wildcards, input, attempt: input.size_mb * 15 + 55000 #+ (attempt - 1) * 200000
    threads: 2
    conda:
        'cp4'
        #'../envs/cellprofiler.yaml'
    #retries: 2
    shell:
        '{params.cellprofiler_executable} -c -r -p {input.pipeline} -i ' + phenotyping_dir + '{wildcards.path}/cellprofiler{wildcards.cycle} -o ' + phenotyping_dir + '{wildcards.path}/cellprofiler{wildcards.cycle}/{wildcards.pipeline}'
        # This command will ignore error from cellprofiler, sometimes necessary if there is a bug:
        #'cellprofiler -c -r -p {input.pipeline} -i ' + phenotyping_dir + '{wildcards.path}/cellprofiler -o ' + phenotyping_dir + '{wildcards.path}/cellprofiler/{wildcards.pipeline} || (test $? = 1 -o $? = 137 && echo \'""\' > {output.cell_file} )'
        #'~/miniconda3/envs/cp4/bin/cellprofiler -c -r -p {input.pipeline} -i ' + phenotyping_dir + '{wildcards.path}/cellprofiler -o ' + phenotyping_dir + '{wildcards.path}/cellprofiler/{wildcards.pipeline}'

rule copy_cellprofiler_output:
    input:
        cell_file = phenotyping_dir + '{path}/cellprofiler{cycle}/{pipeline}/Cells.csv'
    output:
        data = temp(phenotyping_dir + '{path}/cellprofiler{cycle,|cycle\d+}_{pipeline,[^./]+}.csv'),
    run:
        import pandas
        '''
        command = '~/miniconda3/envs/cp4/bin/cellprofiler -c -r -p {pipeline_file} -i {proc}{path}/cellprofiler -o {proc}{path}/cellprofiler/{pipeline}'
        command = command.format(proc=phenotyping_dir, path=wildcards.path, pipeline=wildcards.pipeline, pipeline_file=input.pipeline)
        print (command)
        status = os.system(command)
        code = os.waitstatus_to_exitcode(status)
        debug ("Exit code", code)
        assert code not in (137,), "Cellprofiler failed in an unexpected way"

        table_path = phenotyping_dir + wildcards.path + '/cellprofiler/' + wildcards.pipeline + '/Cells.csv'
        if code == 0:
            table = pandas.read_csv(table_path, index_col=0)
        else:
            table = pandas.DataFrame()
        #os.system('ln -s "{}" "{}"'.format(os.path.relpath(table_path, os.path.dirname(output.data)), output.data))
        '''
        table = pandas.read_csv(input.cell_file, index_col=0)
        table.to_csv(output.data)

        #os.system('touch {}'.format(output.mark))


##################################################
## Custom segmentation for LMNA phenotyping
##################################################

rule run_special_segmentation:
    input:
        segmentation_dir + '{path}/cells.csv',
        get_phenotyping_pt,
        segmentation_dir + '{path}/cells_mask.tif',
        segmentation_dir + '{path}/nuclei_mask.tif',
    output:
        phenotyping_dir + '{path}/puncta_mask.tif',
        phenotyping_dir + '{path}/line_mask.tif',
    resources:
        cuda = 1,
        mem_mb = lambda wildcards, input, attempt: input.size_mb * 2 + 4000,
    run:
        command = '~/miniconda3/envs/ai/bin/python3 segment_lmna.py {} {}'.format(' '.join(input), ' '.join(output))
        code = os.system(command)
        assert code == 0

##################################################
## Merging phenotype tables
##################################################

def get_other_tables(wildcards):
    return [phenotyping_dir + '{path}/' + name + '.csv' for name in wildcards.phenotype_tables.split('.')]

rule merge_tables_phenotype:
    input:
        cell_table = segmentation_dir + '{path}/cells.csv',
        other_tables = get_other_tables,
    output:
        table = temp(phenotyping_dir + '{path}/{phenotype_tables}.cells_phenotype.csv'),
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 2 + 10000
    run:
        import pandas

        cell_table = pandas.read_csv(input.cell_table, index_col=0)
        index = cell_table.index
        #cell_table = cell_table.reset_index(drop=True)
        #table = pandas.concat([cell_table] + [pandas.read_csv(path).iloc[:len(cell_table.index),:] for path in input.other_tables], axis=1)
        table = pandas.concat([pandas.read_csv(path).iloc[:len(cell_table.index),:] for path in input.other_tables], axis=1)
        table = table.set_index(index[:len(table)])
        table.to_csv(output.table)

def get_grid_filenames_pheno(wildcards):
    grid_size = int(wildcards.grid_size)
    numbers = ['{:02}'.format(i) for i in range(grid_size)]
    return expand(phenotyping_dir + '{well}_grid{grid_size}/tile{x}x{y}y/{type}.cells_phenotype.csv', x=numbers, y=numbers, allow_missing=True)

rule merge_grid_pheno_tables:
    input:
        tables = get_grid_filenames_pheno,
        #composite = stitching_dir + '{well}_grid{grid_size}/grid_composite.json',
    output:
        table = temp(phenotyping_dir + '{well}_grid{grid_size,\d+}/{type,[^/]*}.cells_phenotype.csv'),
    resources:
        #mem_mb = lambda wildcards, input: input.size_mb * 50 + 5000
        mem_mb = 5000
    wildcard_constraints:
        #possible_output_dir = '(' + output_dir + ')|',
        type = '[^/]+',
    run:
        def row_func(row):
            i = row['file_index']
            row['pheno_file_path'] = row['file_path']
            row['pheno_tile_index'] = i
            row['pheno_tile_x'] = i // int(wildcards.grid_size)
            row['pheno_tile_y'] = i % int(wildcards.grid_size)

        merge_csv_files(input.tables, output.table, extra_columns=['pheno_file_path', 'pheno_tile_x', 'pheno_tile_y', 'pheno_tile_index'], row_func=row_func)


phenotyping_grid_size = config.get('phenotyping_grid_size', 1)

rule link_merged_grid_phenotype:
    input:
        ((phenotyping_dir + '{well}_grid' + str(phenotyping_grid_size) + '/{type}.cells_phenotype.csv')
                if phenotyping_grid_size != 1 else
                (phenotyping_dir + '{well}/{type}.cells_phenotype.csv')),
    output:
        temp(phenotyping_dir + '{well}_grid/{type,[^/]*}.cells_phenotype.csv'),
    localrule: True
    wildcard_constraints:
        type = '[^/]+',
    shell:
        "cp -l {input[0]} {output[0]}"

ruleorder: link_merged_grid_phenotype > merge_tables_phenotype


