import os
import sys
import glob

##################################################
## Phenotyping
##################################################

def get_phenotyping_pt(wildcards):
    if config['phenotyping']['use_corrected']:
        return stitching_dir + '{path}/corrected_pt.tif'
    else:
        return stitching_dir + '{path}/raw_pt.tif'

rule calc_features:
    input:
        cell_table = segmentation_dir + '{path}/cells.csv',
        cells = segmentation_dir + '{path}/cells_mask.tif',
        nuclei = segmentation_dir + '{path}/nuclei_mask.tif',
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

rule copy_cellprofiler_files:
    input:
        #image = stitching_dir + '{path}/cycle' + phenotype_cycle + '.tif',
        image = get_phenotyping_pt,
        cells = segmentation_dir + '{path}/cells_mask.tif',
        nuclei = segmentation_dir + '{path}/nuclei_mask.tif',
        #puncta = phenotyping_dir + '{path}/puncta_mask.tif',
        #lines = phenotyping_dir + '{path}/line_mask.tif',
    output:
        file_list = phenotyping_dir + '{path}/cellprofiler{cycle}/files.csv',
        images = expand(phenotyping_dir + '{path}/cellprofiler{cycle}/channel{channel}.tif', channel=range(len(config['phenotyping_channels'])), allow_missing=True),
        cells = temp(phenotyping_dir + '{path}/cellprofiler{cycle}/cells.tif'),
        nuclei = temp(phenotyping_dir + '{path}/cellprofiler{cycle}/nuclei.tif'),
        #puncta = phenotyping_dir + '{path}/cellprofiler/puncta.tif',
        #lines = phenotyping_dir + '{path}/cellprofiler/lines.tif',
    wildcard_constraints:
        cycle = '|cycle\d+',
    params:
        cycle = parse_param('cycle', None)
    run:
        import numpy as np
        import tifffile

        image = tifffile.imread(input.image)
        if params.cycle is not None:
            image = image[params.cycle]
        else:
            image = image.reshape(-1, image.shape[2], image.shape[3])

        #remove this
        #bad_shape = tifffile.memmap(input.cells, mode='r').shape

        with open(output.file_list, 'w') as ofile:
            if len(input) > 3:
                ofile.write(','.join(['FileName_CH{}'.format(i) for i in range(len(config['phenotyping_channels']))]) + ',FileName_Cells,FileName_Nuclei,FileName_Puncta,FileName_Line\n')
            else:
                ofile.write(','.join(['FileName_CH{}'.format(i) for i in range(len(config['phenotyping_channels']))]) + ',FileName_Cells,FileName_Nuclei\n')
            for i,path in enumerate(output.images):
                ofile.write(os.path.basename(path) + ',')
                tifffile.imwrite(path, image[i])
                #tifffile.imwrite(path, image[i,:bad_shape[0],:bad_shape[1]])

            for path in input[1:]:
                with tifffile.TiffFile(path) as cells_file:
                    dtype = cells_file.pages[0].dtype
                    assert dtype == np.uint16, (
                        "Segmentation input to cellprofiler must be in uint16 form. "
                        "'{}' has dtype {}, this may be due to there being more than 65535 "
                        "cells in the image. To solve this increase the grid size for phenotyping "
                        " in config.yaml.".format(path, dtype))

            #os.symlink(input.cells, output.cells)
            os.symlink(os.path.relpath(input.cells, os.path.dirname(output.cells)), output.cells)
            ofile.write(os.path.basename(output.cells))
            os.symlink(os.path.relpath(input.nuclei, os.path.dirname(output.nuclei)), output.nuclei)
            ofile.write(',' + os.path.basename(output.nuclei))

            if len(input) > 3:
                os.symlink(os.path.relpath(input.puncta, os.path.dirname(output.puncta)), output.puncta)
                ofile.write(',' + os.path.basename(output.puncta))
                os.symlink(os.path.relpath(input.lines, os.path.dirname(output.lines)), output.lines)
                ofile.write(',' + os.path.basename(output.lines))

            ofile.write('\n')

def find_pipeline(wildcards):
    pipeline = glob.glob('*{}.cppipe'.format(wildcards.pipeline))
    if len(pipeline) != 1:
        return wildcards.pipeline + '.cppipe'
    return pipeline[0]

rule run_cellprofiler:
    input:
        file_list = phenotyping_dir + '{path}/cellprofiler{cycle}/files.csv',
        images = expand(phenotyping_dir + '{path}/cellprofiler{cycle}/channel{channel}.tif', channel=range(len(config['phenotyping_channels'])), allow_missing=True),
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
        mem_mb = lambda wildcards, input, attempt: input.size_mb * 150 + 55000 #+ (attempt - 1) * 200000
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
        table = phenotyping_dir + '{path}/{phenotype_tables}.cells_phenotype.csv',
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
        table = phenotyping_dir + '{well}_grid{grid_size,\d+}/{type,[^/]*}.cells_phenotype.csv',
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
        phenotyping_dir + '{well}_grid/{type,[^/]*}.cells_phenotype.csv',
    localrule: True
    wildcard_constraints:
        type = '[^/]+',
    shell:
        "cp -l {input[0]} {output[0]}"

ruleorder: link_merged_grid_phenotype > merge_tables_phenotype


