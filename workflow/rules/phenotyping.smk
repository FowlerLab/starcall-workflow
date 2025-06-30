import os
import sys
import glob

##################################################
## Phenotyping
##################################################

rule calc_features:
    input:
        cell_table = phenotyping_input_dir + '{prefix}/cells.csv',
        cells = phenotyping_input_dir + '{prefix}/cells_mask.tif',
        nuclei = phenotyping_input_dir + '{prefix}/nuclei_mask.tif',
        image = phenotyping_input_dir + '{prefix}/corrected_pt.tif',
    output:
        features = phenotyping_output_dir + '{prefix}/features.csv'
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

        for index, cell in starcall.utils.simple_progress(list(cell_table.iterrows())):
            x1, y1 = cell.bbox_x1 * phenotype_scale, cell.bbox_y1 * phenotype_scale
            x2, y2 = cell.bbox_x2 * phenotype_scale, cell.bbox_y2 * phenotype_scale
            cell_mask = cells[x1:x2,y1:y2] == index
            nucleus_mask = nuclei[x1:x2,y1:y2] == index
            image_section = image[:,x1:x2,y1:y2]
            props = all_props[index]

            for prop in props:
                if type(props[prop]) in (float, int):
                    features.setdefault(prop, []).append(props[prop])

            for maskname, mask in [('cell', cell_mask), ('nucleus', nucleus_mask), ('cytoplasm', cell_mask & ~nucleus_mask)]:
                for channel in range(phenotype_channels):
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
        #image = phenotyping_input_dir + '{prefix}/cycle' + phenotype_cycle + '.tif',
        image = phenotyping_input_dir + '{prefix}/corrected_pt.tif',
        cells = phenotyping_input_dir + '{prefix}/cells_mask.tif',
        nuclei = phenotyping_input_dir + '{prefix}/nuclei_mask.tif',
        #puncta = phenotyping_dir + '{prefix}/puncta_mask.tif',
        #lines = phenotyping_dir + '{prefix}/line_mask.tif',
    output:
        file_list = phenotyping_dir + '{prefix}/cellprofiler/files.csv',
        images = expand(phenotyping_dir + '{prefix}/cellprofiler/channel{channel}.tif', channel=range(phenotype_channels), allow_missing=True),
        cells = phenotyping_dir + '{prefix}/cellprofiler/cells.tif',
        nuclei = phenotyping_dir + '{prefix}/cellprofiler/nuclei.tif',
        #puncta = phenotyping_dir + '{prefix}/cellprofiler/puncta.tif',
        #lines = phenotyping_dir + '{prefix}/cellprofiler/lines.tif',
    run:
        import numpy as np
        import tifffile

        image = tifffile.imread(input.image)
        image = image.reshape(-1, image.shape[2], image.shape[3])

        #remove this
        bad_shape = tifffile.memmap(input.cells, mode='r').shape

        with open(output.file_list, 'w') as ofile:
            if len(input) > 3:
                ofile.write(','.join(['FileName_CH{}'.format(i) for i in range(phenotype_channels)]) + ',FileName_Cells,FileName_Nuclei,FileName_Puncta,FileName_Line\n')
            else:
                ofile.write(','.join(['FileName_CH{}'.format(i) for i in range(phenotype_channels)]) + ',FileName_Cells,FileName_Nuclei\n')
            for i,path in enumerate(output.images):
                ofile.write(os.path.basename(path) + ',')
                #tifffile.imwrite(path, image[i])
                tifffile.imwrite(path, image[i,:bad_shape[0],:bad_shape[1]])

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
        file_list = phenotyping_dir + '{prefix}/cellprofiler/files.csv',
        images = expand(phenotyping_dir + '{prefix}/cellprofiler/channel{channel}.tif', channel=range(phenotype_channels), allow_missing=True),
        cells = phenotyping_dir + '{prefix}/cellprofiler/cells.tif',
        nuclei = phenotyping_dir + '{prefix}/cellprofiler/nuclei.tif',
        #puncta = phenotyping_dir + '{prefix}/cellprofiler/puncta.tif',
        #lines = phenotyping_dir + '{prefix}/cellprofiler/lines.tif',
        pipeline = find_pipeline,
        #pipeline = '{pipeline}.cppipe',
    output:
        #data = phenotyping_output_dir + '{prefix}/cellprofiler_{pipeline,[^./]+}.csv',
        #mark = phenotyping_dir + '{prefix}/cellprofiler/{pipeline,[^./]+}/mark',
        cell_file = phenotyping_dir + '{prefix}/cellprofiler/{pipeline}/Cells.csv'
    resources:
        mem_mb = lambda wildcards, input, attempt: input.size_mb * ([50, 100, 500][attempt-1]) + 10000 #+ (attempt - 1) * 200000
    threads: 2
    conda:
        'cp4'
        #'../envs/cellprofiler.yaml'
    #retries: 2
    shell:
        'cellprofiler -c -r -p {input.pipeline} -i ' + phenotyping_dir + '{wildcards.prefix}/cellprofiler -o ' + phenotyping_dir + '{wildcards.prefix}/cellprofiler/{wildcards.pipeline} || test $? = 0 -o test $? = 1 && echo '""' > {output.cell_file}'
        #'~/miniconda3/envs/cp4/bin/cellprofiler -c -r -p {input.pipeline} -i ' + phenotyping_dir + '{wildcards.prefix}/cellprofiler -o ' + phenotyping_dir + '{wildcards.prefix}/cellprofiler/{wildcards.pipeline}'

rule copy_cellprofiler_output:
    input:
        cell_file = phenotyping_dir + '{prefix}/cellprofiler/{pipeline}/Cells.csv'
    output:
        data = phenotyping_output_dir + '{prefix}/cellprofiler_{pipeline,[^./]+}.csv',
    run:
        import pandas
        '''
        command = '~/miniconda3/envs/cp4/bin/cellprofiler -c -r -p {pipeline_file} -i {proc}{prefix}/cellprofiler -o {proc}{prefix}/cellprofiler/{pipeline}'
        command = command.format(proc=phenotyping_dir, prefix=wildcards.prefix, pipeline=wildcards.pipeline, pipeline_file=input.pipeline)
        print (command)
        status = os.system(command)
        code = os.waitstatus_to_exitcode(status)
        debug ("Exit code", code)
        assert code not in (137,), "Cellprofiler failed in an unexpected way"

        table_path = phenotyping_dir + wildcards.prefix + '/cellprofiler/' + wildcards.pipeline + '/Cells.csv'
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
        phenotyping_input_dir + '{prefix}/cells.csv',
        phenotyping_input_dir + '{prefix}/corrected_pt.tif',
        phenotyping_input_dir + '{prefix}/cells_mask.tif',
        phenotyping_input_dir + '{prefix}/nuclei_mask.tif',
    output:
        phenotyping_dir + '{prefix}/puncta_mask.tif',
        phenotyping_dir + '{prefix}/line_mask.tif',
    resources:
        cuda = 1,
        mem_mb = lambda wildcards, input, attempt: input.size_mb * 2 + 4000,
    run:
        command = '~/miniconda3/envs/ai/bin/python3 segment_lmna.py {} {}'.format(' '.join(input), ' '.join(output))
        code = os.system(command)
        assert code == 0

