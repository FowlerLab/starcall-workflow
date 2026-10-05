import os
import glob
import re
import pandas as pd
from collections import OrderedDict


def get_aux_data_correction_summary(wildcards):
    wildcards.path = wildcards.well + '_grid' + wildcards.grid_size
    return get_aux_data(wildcards)


def get_aux_data(wildcards, path=None):
    path = wildcards.path if path is None else path

    #if path != '': path = path + '.'

    files = []
    for base_dir in (sequencing_dir, input_dir):
        pattern = base_dir + '{path}/{segmentation_type}.auxdata/*.csv'.format(path=path, segmentation_type=wildcards.segmentation_type)
        files.extend(sorted(glob.glob(pattern)))
        pattern = base_dir + '{path}/auxdata/*.csv'.format(path=path, segmentation_type=wildcards.segmentation_type)
        files.extend(sorted(glob.glob(pattern)))

    #if re.fullmatch('(tile.+)|(well.+)|(cycle.+)', os.path.basename(path)):
    if path.count('_grid'):
        files.extend(get_aux_data(wildcards, path=path.split('_grid')[0]))
    elif path:
        files.extend(get_aux_data(wildcards, path=os.path.dirname(path)))

    for i in range(len(files)):
        files[i] = files[i].replace('//', '/')

    return files



def find_dots_mem_mb(wildcards, input, threads):
    if config['dotdetection'].get('backend', 'original') == 'original':
        return size_mb(input) * 10 + 15000
    # one float32 copy of the sequencing channels, plus per thread bands and tiles
    return size_mb(input) * 4 + 10000 + 500 * threads


rule find_dots:
    """ Detect amplicon colonies in the sequencing images.
    This is a crutial step in sequencing the barcodes expressed in cells, and has
    a couple steps:
        First, background and cell debris is filtered with a difference of gaussian filter,
        Then, differing intensities of sequencing channels are corrected by z-score normalizing
            each channel and cycle
        Then dots are highlighted by subtracting the second maximal channel from all channels on
            a per-pixel basis, and all cycles are combined by taking the standard deviation across
            cycles, again on a per-pixel basis.
        This greyscale image is given to the laplacian of gaussian blob detection algorith
            (https://scikit-image.org/docs/stable/api/skimage.feature.html#skimage.feature.blob_log)
            with the parameters min, max, num, specifying the gaussian sigmas to search over
        This results in an x,y position for each detected colony. We extract read values from this position,
            and save these values to the output csv file.

    Params:
        min, max, mean: The range of gaussian sigmas to search for colonies. 1-3 should capture most colonies.

    Config:
        dotdetection.backend: original, cpu or gpu, which version of detect_dots to run (see above)
        dotdetection.threads: threads to use, only the cpu and gpu versions use more than one

    Output:
        The output is a csv file containing the columns:
            position_x, position_y: The pixel position of the colony
            values_cycle00_G, values_cycle00_T, ...:
                The values extracted from the sequencing images at the colony position,
                for each cycle and channel.
    """
    input:
        sequencing_dir + '{path}/raw.tif',
    output:
        sequencing_dir + '{path}/bases{min}{max}{num}.csv',
        #sequencing_dir + '{path}/dot_filter.tif',
        **({
            'diagnostics_summary': sequencing_dir + '{path}/dotdiagnostics{min}{max}{num}.summary.csv',
            'diagnostics_images': sequencing_dir + '{path}/dotdiagnostics{min}{max}{num}.images.npz',
        } if config['dotdetection'].get('diagnostics', False) else {}),
    params:
        min_sigma = parse_param('min', config['dotdetection']['min_sigma']),
        max_sigma = parse_param('max', config['dotdetection']['max_sigma']),
        num_sigma = parse_param('num', config['dotdetection']['num_sigma']),
        backend = config['dotdetection'].get('backend', 'original'),
        #settings for the sample to build for the quality report site
        diagnostics_stride = config['dotdetection'].get('pixel_sample_stride', 8),
        diagnostics_crop_size = config['dotdetection'].get('crop_sample_size', 256),
        diagnostics_dots = config['dotdetection'].get('max_dots_for_sample', 50000),

    wildcard_constraints:
        min = '|_min\d+(.\d+)?',
        max = '|_max\d+(.\d+)?',
        num = '|_num\d+(.\d+)?',
        
    resources:
        mem_mb = find_dots_mem_mb,
        cuda = 1 if config['dotdetection'].get('backend', 'original') == 'gpu' else 0,
    threads: config['dotdetection'].get('threads', 4 if config['dotdetection'].get('backend', 'original') == 'original' else 16)
    run:
        import numpy as np
        import tifffile
        import starcall.dotdetection
        import starcall.correction
        import skimage.morphology

        full_well = tifffile.memmap(input[0], mode='r')
        if params.backend == 'original':
            image = full_well[...,sequencing_channels_slice,:,:].astype(np.float32, copy=True)
            del full_well
            detect_dots, kwargs = starcall.dotdetection.detect_dots, {}
        elif params.backend in ('cpu', 'gpu'):
            # these read bands from the memmap as needed, the input is never modified
            image = full_well[...,sequencing_channels_slice,:,:]
            detect_dots = getattr(starcall.dotdetection, 'detect_dots_' + params.backend)
            kwargs = dict(threads = threads)
            if 'diagnostics_summary' in output.keys():
                kwargs['diagnostics'] = output.diagnostics_summary[:-len('.summary.csv')]
                kwargs['diagnostics_stride'] = params.diagnostics_stride
                kwargs['diagnostics_crop_size'] =  params.diagnostics_crop_size
                kwargs['diagnostics_max_sample_dots'] = params.diagnostics_dots
        else:
            raise ValueError("dotdetection.backend should be original, cpu or gpu, not '{}'".format(params.backend))

        if not np.any(image):
            reads = pandas.DataFrame()
            if 'diagnostics_summary' in output.keys():
                # empty placeholders, skipped when plotting
                import starcall.dotdetection_qc
                pandas.DataFrame(columns=starcall.dotdetection_qc.SUMMARY_COLUMNS).to_csv(output.diagnostics_summary, index=False)
                np.savez_compressed(output.diagnostics_images)
        else:
            debug('keeping z-scored intensities...')
            debug('running dot detection with backend {} and {} threads'.format(params.backend, threads))
            reads = detect_dots(
                image,
                min_sigma = params.min_sigma,
                max_sigma = params.max_sigma,
                num_sigma = params.num_sigma,
                copy = False,
                channels = sequencing_channels_order,
                debug = debug,
                **kwargs,
            )

        reads.to_csv(output[0])


def get_cycle_str(i):
    if i >= 10:
        return str(i)
    else:
        return "0" + str(i)

rule attach_quality_information:
    """ Attach the PhredQ like score to the reads, along with  """
     input:
        sequencing_dir + '{path}/bases{params}.csv',
    output:
        sequencing_dir + '{path}/quality_bases{params}.csv',
    wildcard_constraints:
        params = params_regex('min', 'max', 'num'),
    resources:
        mem_mb = lambda wildcards, input: size_mb(input) * 3 + 15000
    run:
        import pandas as pd
        from starcall.qc import get_softmax_df, calculate_peaks, get_chastity_df, get_log_margin_df, get_purity_df, get_dominance_signed_df, get_fac_delta_of_top_pos
        import numpy as np
        import starcall.reads #does this change what is called?

        reads = pd.read_csv(input[0], index_col=0)
        orig_rows = reads.shape[0]

        #filter reads table to remove nan rows (any dots which are too close to the edge are forced to NaN)
        reads = reads[~reads.isna().any(axis=1)].copy()

        debug ('removed ', orig_rows - reads.shape[0], ' which are NaNs due to proximity to edges')

        
        value_cols = [c for c in reads.columns if c.startswith('values_cycle')]

        #attach phredq score and intensity change info
        quality_scores = get_softmax_df(reads) #use z scored values 
        peaks = calculate_peaks(reads) #use original values
        #testing other metrics 
        chastity_scores = get_chastity_df(reads)
        for i in range(0, quality_scores.shape[-1]):
            reads['chastity_cycle' + get_cycle_str(i)] = chastity_scores[:, i]
        reads['mean_chastity'] = np.mean(chastity_scores, axis=1)
        reads['min_chastity'] = np.min(chastity_scores, axis=1)
        reads['max_seq'] = reads.reads.sequences 

        reads.to_csv(output[0])


rule attach_cell_ids:
    """ Attaches the ID of each cell the dots are in to said dots

    Output: The output is a csv file containing the columns:
            position_x, position_y: The pixel position of the colony
            values_cycle00_G, values_cycle00_T, ...:
                The values extracted from the sequencing images at the colony position,
                for each cycle and channel.
            cell: The cell that each read is contained in, 0 if not in a cell.
    """
    input:
        bases = sequencing_dir + '{path}/quality_bases{params}.csv',
        cells = segmentation_dir + '{path}/{segmentation_type}_mask_downscaled.tif',
        #cells_table = segmentation_dir + '{path}/{segmentation_type}.csv',
    output:
        table = sequencing_dir + '{path}/{segmentation_type}_quality{params}.csv',
    wildcard_constraints:
        params = params_regex('min', 'max', 'num'),
    resources:
        mem_mb = lambda wildcards, input: 5000 +  size_mb(input) * 10
    run:
        import tifffile
        import numpy as np
        import pandas
        import csv
        import matplotlib.pyplot as plt
        import starcall.reads

        cells = tifffile.imread(input.cells)
        table = pandas.read_csv(input.bases, index_col=0)
        xposes, yposes = np.round(table.reads.positions.T).astype(int)
        table['cell'] = cells[xposes,yposes]
        table.to_csv(output.table)



def get_orig_method_files(wildcards):
    grid_size = int(wildcards.grid_size)
    numbers = ['{:02}'.format(i) for i in range(grid_size)]
    return expand(sequencing_dir + '{well}_grid{grid_size}/tile{x}x{y}y/cells_quality.csv', x=numbers, y=numbers, allow_missing=True)

#TODO -  move this to QC and make this into a table of qcs per subtile
rule make_orig_method_comparison_table:
    input:
        orig_tables = get_orig_method_files,
        library = get_aux_data_correction_summary,
    output:
        summ_table = sequencing_dir + '{well}_grid{grid_size}/{segmentation_type}_summary_orig_approach.csv',
    resources:
        mem_mb = lambda wildcards, input: 5000 +  size_mb(input) * 2
    run:    
        import pandas as pd 
        import re

        #setup exact match table
        barcodes_table = pd.read_csv(input.library[0])
        #after barcode table corrections, all barcodes for matching should be length 12 
        #and the first columns should contain the sequences to be matched with
        barcode_table_cols = list(barcodes_table.columns)
        dummy_barcodes2 = barcodes_table[barcode_table_cols[:1]].rename(columns = {barcode_table_cols[0]: 'orig_barcode_match'})

        summary_rows = []
        for table_path in input.orig_tables:
            tile_match = re.search(r'/(tile\d+x\d+y)/', table_path)
            tile = tile_match.group(1) if tile_match else table_path

            table = pd.read_csv(table_path, index_col = 0)
            #filter to only those in cells 
            table = table[table.cell != 0].copy()
            #merge exact matches for barcodes
            table = table.merge(dummy_barcodes2, left_on = 'max_seq', right_on = 'orig_barcode_match', how = 'left')

            orig_matched = ~table['orig_barcode_match'].isna()
            

            summary_rows.append({
                'tile': tile,
                'n_reads': len(table),
                'matched_reads': int(orig_matched.sum()),
                'mean_phred_before': table['mean_chastity'].mean(),
                'median_phred_before': table['min_chastity'].median(),
                'mean_min_phred_before': table['min_phred'].mean(),
                'median_min_phred_before': table['min_phred'].median(),
               
            })

        summary_table = pd.DataFrame(summary_rows)
        summary_table.to_csv(output.summ_table)


rule make_orig_method_metric_performance_table:
    input:
        orig_tables = get_orig_method_files, #cells quality files have all the metrics attached
        library = get_aux_data_correction_summary,
    output:
        summ_table = sequencing_dir + '{well}_grid{grid_size}/{segmentation_type}_summary_quality_metrics.csv',
    resources:
        mem_mb = lambda wildcards, input: 5000 +  size_mb(input) * 2
    run:    
        import pandas as pd 
        import re
        from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve

        #setup exact match table
        barcodes_table = pd.read_csv(input.library[0])
        #after barcode table corrections, all barcodes for matching should be length 12 
        #and the first columns should contain the sequences to be matched with
        barcode_table_cols = list(barcodes_table.columns)
        dummy_barcodes2 = barcodes_table[barcode_table_cols[:1]].rename(columns = {barcode_table_cols[0]: 'orig_barcode_match'})

        metric_cols = [
                'mean_phred', 'min_phred', 'sum_std_intensities',
                'mean_chastity', 'min_chastity',
                'mean_log_margin', 'min_log_margin',
                'mean_purity', 'min_purity',
                'mean_dominance_signed', 'min_dominance_signed',
                'min_frac_delta_of_top', 'mean_frac_delta_of_top']

        summary_rows = []
        for table_path in input.orig_tables:
            tile_match = re.search(r'/(tile\d+x\d+y)/', table_path)
            tile = tile_match.group(1) if tile_match else table_path

            table = pd.read_csv(table_path, index_col = 0)
            #filter to only those in cells 
            table = table[table.cell != 0].copy()
            #merge exact matches for barcodes
            table = table.merge(dummy_barcodes2, left_on = 'max_seq', right_on = 'orig_barcode_match', how = 'left')

            table['has_match'] = ~table.orig_barcode_match.isna()
            y_all = table['has_match'].to_numpy()
            
            for col in metric_cols:
                x_all = table[col].to_numpy()
                valid = ~pd.isna(x_all)
                x, y = x_all[valid], y_all[valid]

                auc = roc_auc_score(y, x)
                ap = average_precision_score(y, x)

                summary_rows.append({
                    'metric': col,
                    'tile': tile, 
                    'auroc': auc,                 
                    'avg_precision': ap,  
                    'n_match': int(y.sum()),
                    'n_nomatch': int((~y).sum()),         

                })

        summary_table = pd.DataFrame(summary_rows)
        summary_table.to_csv(output.summ_table)


ruleorder: segment_cells > segment_cells_bases





rule annotate_dots:
    """ Marks each dot that was detected with a cross. The annotations are included in a new
    channel added to the image
    """
    input:
        image = sequencing_dir + '{path}/raw.tif',
        bases = sequencing_dir + '{path}/bases{params}.csv',
    output:
        qc_dir + '{path}/annotated{params}.tif',
    wildcard_constraints:
        params = params_regex('min', 'max', 'num'),
    resources:
        mem_mb = 25000
    run:
        import numpy as np
        import starcall.utils
        import tifffile
        import skimage.draw
        import starcall.reads
        import pandas

        table = pandas.read_csv(input.bases, index_col=0)
        image = tifffile.memmap(input[0], mode='r')[0]
        marked_image = starcall.utils.mark_dots(image, table.reads.positions.astype(int))
        tifffile.imwrite(output[0], marked_image)


rule merge_final_tables:
    """ Combine the read table with any additional tables that have information linked to specific barcodes.
    This could include variants, gene kos, or other perturbations depending on the experiment.

    Additional tables are searched for in the following paths:
        sequencing/{path}/{segmentation_type}.auxdata/
        sequencing/{path}/auxdata/
        input/{path}/{segmentation_type}.auxdata/
        input/{path}/auxdata/
        input/auxdata/
    All tables found should have a barcode column as the first column, containing barcodes that will be matched
    to reads. If multiple reads are required for a match, they should be separated with a dash, for example the
    barcode 'GTAC-AATG' will only be matched to a cell with both 'GTAC' and 'AATG' as reads.

    Read matching is also performed for cells that don't have a perfect match to any barcode, in which case the search
    is expanded to find the barcode that have the minimum edit distance to one of the cells reads. If there are multiple
    such barcodes, no match is able to be made.

    Cells that have been matched to a barcode will have the entire row corresponding to that barcode concatenated
    to the end of the row. Cells that were not matched will not have values for all of these rows, which results in NA
    for all values. This can be prevented by removing all cells that were not able to be matched, enabled in config.yaml
    with the option sequencing.remove_unmatched_cells
    """
    input:
        cell_table = sequencing_dir + '{path}/{segmentation_type}_reads_partial{params}.csv',
        aux_data = get_aux_data,
    output:
        full_table = sequencing_dir + '{path}/{segmentation_type}_reads_old{params}.csv',
    params:
        remove_unmatched = config['sequencing'].get('remove_unmatched_cells', False),
    wildcard_constraints:
        params = params_regex('min', 'max', 'num', 'norm', 'posweight', 'valweight', 'seqweight', 'thresh', 'linkage', 'maxreads'),
    resources:
        mem_mb = lambda wildcards, input: 5000 +  size_mb(input) * 250
    run:
        import pandas
        import numpy as np
        import starcall.sequencing
        import warnings

        debug ('begin ')
        cell_table = pandas.read_csv(input.cell_table, index_col=0)
        debug (cell_table)
        cell_cols = cell_table.columns.copy()

        def join_barcode(cell_table, aux_table):
            if 'read_0' not in cell_table.columns: 
                return cell_table #returning the empty table if no reads are found for this tile - should not affect the merging at the end
            barcodes = []
            for tmp_read in cell_table['read_0']:
                if type(tmp_read) == str: break
            num_cycles = len(tmp_read)
            #num_cycles = len(cell_table['read_0'].iloc[0])
            debug ('num_cycles', num_cycles)
            for barc in aux_table.index:
                new_barcodes = barc.split('-')
                for new_barc in new_barcodes:
                    if len(new_barc) != num_cycles:
                        warnings.warn('Barcode not the right length {} {}'.format(len(new_barc), num_cycles))
                barcodes.append([subbarc[:num_cycles] for subbarc in new_barcodes])

            barcodes = np.array(barcodes)
            lengths = np.array(list(map(len, barcodes.flat)))
            debug (lengths)
            debug (lengths.min(), lengths.mean(), lengths.max())

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
            reads = reads.astype('U' + str(num_cycles))
            barcodes = barcodes.astype('U' + str(num_cycles))
            counts[np.isnan(counts)] = 0
            counts = counts.astype(int)
            debug (reads.dtype)
            debug (barcodes.dtype)

            debug ('starting matching')
            indices, read_indices, edit_distances = starcall.sequencing.match_barcodes(reads, counts, barcodes, n_neighbors=2, max_edit_distance=99999, return_distances=True, debug=True, progress=True)
            #edit_distances, indices = library.nearest(reads, counts)
            debug ('  done')
            debug (np.sum(indices[:,0] != -1) / len(indices))

            multiple_matches = edit_distances[:,0] == edit_distances[:,1]
            debug (multiple_matches.shape)

            if params.remove_unmatched:
                new_rows = [aux_table.iloc[i,:] for i in indices[:,0] if i != -1]
                new_index = [cell_table.index[i] for i in indices[:,0] if i != -1]
            else:
                new_rows = [pandas.Series(dtype=object) if i == -1 else aux_table.iloc[i,:] for i in indices[:,0]]
                new_index = cell_table.index

            #new_rows = [pandas.Series(dtype=object) if (i == -1 and not multiple) else aux_table.iloc[i,:] for i, multiple in zip(indices[:,0], multiple_matches)]
            debug ('making new table', len(new_rows))
            new_table = pandas.DataFrame(new_rows, index=new_index)
            debug ('  done')
            new_table['editDistance'] = edit_distances[:,0]
            new_table['edit_distance'] = edit_distances[:,0]
            new_table['matched_barcode_index'] = indices[:,0]
            new_table['edit_distance_2'] = edit_distances[:,1]
            new_table['matched_barcode_index_2'] = indices[:,1]
            debug ('adding rows to table', read_indices.shape[2])
            for i in range(read_indices.shape[2]):
                debug (read_indices[:,0,i].shape)
                new_table['matched_read_index_{}'.format(i)] = read_indices[:,0,i]
                debug ( reads[list(range(len(reads))),read_indices[:,0,i]].shape)

                matched_reads = reads[list(range(len(reads))),read_indices[:,0,i]]
                matched_reads[indices[:,0]==-1] = ''
                new_table['matched_read_{}'.format(i)] = matched_reads

                matched_barcodes = barcodes[indices[:,0]][:,i]
                matched_barcodes[indices[:,0]==-1] = ''
                new_table['matched_barcode_{}'.format(i)] = matched_barcodes
            debug ('added rows to table')
            result = cell_table.join(new_table)
            debug ('joined table')
            return result

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

        debug ('saving to csv')
        cell_table.to_csv(output.full_table)
        debug ('  done')


##################################################
## Merging read tables
##################################################

def get_grid_filenames_seq(wildcards):
    grid_size = int(wildcards.grid_size)
    numbers = ['{:02}'.format(i) for i in range(grid_size)]
    return expand(sequencing_dir + '{well}_grid{grid_size}/tile{x}x{y}y/{segmentation_type}{qc}_reads.csv', x=numbers, y=numbers, allow_missing=True)

rule merge_grid_read_tables:
    input:
        tables = get_grid_filenames_seq,
        composite = lambda wildcards: [stitching_dir + '{well}_grid{grid_size}/grid_composite.json'] if wildcards.qc != '' else [],
    output:
        table = sequencing_dir + '{well}_grid{grid_size,\d+}/{segmentation_type}{qc,|_qc}_reads.csv',
    resources:
        #mem_mb = lambda wildcards, input: input.size_mb * 50 + 5000
        mem_mb = 10000
    run:
        import constitch

        if len(input.composite) != 0:
            composite = constitch.load(input.composite[0])

        def row_func(row):
            i = row['file_index']
            row['seq_file_path'] = row['file_path']
            row['seq_tile_index'] = i
            row['seq_tile_x'] = i // int(wildcards.grid_size)
            row['seq_tile_y'] = i % int(wildcards.grid_size)

            if 'position_x' in row:
                row['position_x'] += composite.boxes[i].position[0]
                row['position_y'] += composite.boxes[i].position[1]

        merge_csv_files(input.tables, output.table, extra_columns=['seq_file_path', 'seq_tile_index', 'seq_tile_x', 'seq_tile_y'], row_func=row_func)


sequencing_grid_size = config.get('sequencing_grid_size', 1)

rule link_merged_grid_reads:
    input:
        ((sequencing_dir + '{well}_grid' + str(sequencing_grid_size) + '/{segmentation_type}_reads.csv')
                if sequencing_grid_size != 1 else
                (sequencing_dir + '{well}/{segmentation_type}_reads.csv')),
    output:
        sequencing_dir + '{well}_grid/{segmentation_type}_reads.csv',
    localrule: True
    wildcard_constraints:
        path_nogrid = '((?!_grid)[^.])*',
    shell:
        "cp -l {input[0]} {output[0]}"

ruleorder: link_merged_grid_reads > merge_final_tables
