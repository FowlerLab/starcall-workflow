""" Rules for matching called sequences against the barcode library.

Collects every distinct sequence across all tiles in a well once and aligns 
them with the barcode table.  Then uses the aligned table to map the closest barcodes per 
read before gathering them into per cell tables with reads ordered by chastity
"""

ALIGN_BATCH_SIZE = 5000

def get_grid_filenames_zscore_with_cells(wildcards):
    grid_size = int(wildcards.grid_size)
    numbers = ['{:02}'.format(i) for i in range(grid_size)]
    return expand(sequencing_dir + '{well}_grid{grid_size}/tile{x}x{y}y/{segmentation_type}_quality{params}.csv', x=numbers, y=numbers, allow_missing=True)

rule collect_unique_sequences_zscored:
    """ Gathers every distinct corrected_max_seq across all tiles in a well/method/approach,
    so the (much slower) barcode alignment step only ever runs once per distinct sequence
    for the whole well, not once per tile.
    """
    input:
        corrected_tables = get_grid_filenames_zscore_with_cells,
    output:
        sequences = sequencing_dir + '{well}_grid{grid_size}/{segmentation_type}_zscore_unique_sequences{params}.csv',
    wildcard_constraints:
        params = params_regex('min', 'max', 'num',),
        grid_size = '\d+',
    resources:
        mem_mb = lambda wildcards, input: 5000 + size_mb(input) * 2
    run:
        import pandas as pd
        import starcall.reads

        all_seqs = []
        for path in input.corrected_tables:
            debug('reading sequences from...', path)
            table = pd.read_csv(path, index_col = 0)
            #max_seq is filled in the quality rule
            all_seqs.append(table['max_seq'])

        all_seqs = pd.concat(all_seqs, axis=0, ignore_index=True)
        unique_seqs = pd.unique(all_seqs.dropna())
        debug (len(all_seqs), 'reads,', len(unique_seqs), 'unique sequences across the well')

        pd.DataFrame({'sequence': unique_seqs}).to_csv(output.sequences, index=False)


rule collect_unique_sequences_zscored_no_grid:
    """ Full well equivalnet of the tiled version - for runs requesting full well values without tiling
    """
    input:
        fullwell_table = sequencing_dir + '{well}/{segmentation_type}_quality{params}.csv'
    output:
        sequences = sequencing_dir + '{well}/{segmentation_type}_zscore_unique_sequences{params}.csv',
    wildcard_constraints:
        params = params_regex('min', 'max', 'num',),
    resources:
        mem_mb = lambda wildcards, input: 5000 + size_mb(input) * 2
    run:
        import pandas as pd
        import starcall.reads

        all_seqs = []
      
        table = pd.read_csv(input.fullwell_table, index_col = 0)
        unique_seqs = table['max_seq'].dropna().unique()
        debug (len(table), 'reads,', len(unique_seqs), 'unique sequences across the well')
        pd.DataFrame({'sequence': unique_seqs}).to_csv(output.sequences, index=False)



rule align_unique_well_sequences_to_barcodes_zscored:
    """ For every unique sequence, finds the barcode(s) achieving that sequence's own true
    minimum Hamming distance to the library (not a fixed cutoff or top-K). Matching is 
    exact-brute-force (one-hot encode + full pairwise Manhattan distance equals Hamming distance here)
    """
    input:
        sequences = sequencing_dir + '{path}/{segmentation_type}_zscore_unique_sequences{params}.csv',
        library = get_aux_data,
    output:
        matches = sequencing_dir + '{path}/{segmentation_type}_zscore_barcode_matches{params}.csv',
    wildcard_constraints:
        params = params_regex('min', 'max', 'num', ),
    params:
        batch_size = ALIGN_BATCH_SIZE,
        num_cycles = config['cycles'],
    threads: config['sequencing'].get('number_cores_for_alignment', 8),
    resources:
        mem_mb = lambda wildcards, input, threads: 5000 + (count_lines(input.library[0]) * ALIGN_BATCH_SIZE * 12 * threads) // 1_000_000
    run:
        import os

        #pin each worker to a single thread -- otherwise sklearn/BLAS spawns its own
        #internal threads per worker process
        os.environ['OMP_NUM_THREADS'] = '1'
        os.environ['OPENBLAS_NUM_THREADS'] = '1'
        os.environ['MKL_NUM_THREADS'] = '1'

        import pandas as pd
        import numpy as np
        from starcall.sequencing import init_barcode_match_worker, match_barcode_batch
        from concurrent.futures import ProcessPoolExecutor
        import os
        import warnings

        #adapted from old sequencing merge_final_tables barcode section
        aux_table = pd.read_csv(input.library[0])
        aux_table.set_index(aux_table.columns[0], inplace = True)
        num_cycles = len(params.num_cycles)
        barcodes = []
        barcodes_full = []
        for barc in aux_table.index:
            new_barcodes = barc.split('-')
            for new_barc in new_barcodes:
                if len(new_barc) != num_cycles:
                    warnings.warn('Barcode not the right length {} {}'.format(len(new_barc), num_cycles))
            barcodes.append([subbarc[:num_cycles] for subbarc in new_barcodes])
            barcodes_full.append(new_barcodes)
        barcode_seqs = np.array(barcodes)

        #sequences_to_vector needs a fixed-width unicode dtype -- this truncated version is
        #used only for the actual distance matching (reads are num_cycles bases long)
        barcode_seqs = np.asarray(barcode_seqs, dtype=f'<U{num_cycles}')
        #untruncated labels (may be longer than num_cycles) used only to build the
        #human-readable barcode_matches output column, so long barcodes aren't clipped
        max_barcode_len = max(len(subbarc) for barc in barcodes_full for subbarc in barc)
        barcode_seqs_full = np.asarray(barcodes_full, dtype=f'<U{max_barcode_len}')
        read_seqs = pd.read_csv(input.sequences)['sequence'].to_list()
        batch_size = params.batch_size

        #check if a version of the alignment table exists already (if it does, add the new sequences only to it
        #this is just to make it easier to prototype sequencing approaches - we can remove this later 
        if os.path.exists(output.matches):
            already_matched = pd.read_csv(output.matches)
            already_matched = already_matched.sequence.to_list()
            debug ('aligmnet table already exists...initial reads to align: ', len(read_seqs))
            read_seqs = set(read_seqs).difference(set(already_matched))
            debug ('reads not already aligned: ', len(read_seqs))
            read_seqs = np.asarray(list(read_seqs), dtype=f'<U{num_cycles}')
            batches = [(start, read_seqs[start:start + batch_size]) for start in range(0, len(read_seqs), batch_size)]

            min_dist = np.empty(len(read_seqs), dtype=int)
            matched_barcodes = [None] * len(read_seqs)

            with ProcessPoolExecutor(max_workers=threads, initializer=init_barcode_match_worker, initargs=(barcode_seqs, barcode_seqs_full)) as executor:
                for start, batch_min, batch_matches in executor.map(match_barcode_batch, batches):
                    min_dist[start:start + len(batch_min)] = batch_min
                    for i, m in enumerate(batch_matches):
                        matched_barcodes[start + i] = m

            result = pd.DataFrame({
                'sequence': read_seqs,
                'min_hamming_distance': min_dist,
                'barcode_matches': matched_barcodes,
            })
            result.to_csv(output.matches, index=False)
        else:
            debug ('aligning', len(read_seqs), 'unique sequences against', len(barcode_seqs), 'barcodes using', threads, 'workers')
           
            read_seqs = np.asarray(list(read_seqs), dtype=f'<U{num_cycles}')
            batches = [(start, read_seqs[start:start + batch_size]) for start in range(0, len(read_seqs), batch_size)]
            min_dist = np.empty(len(read_seqs), dtype=int)
            matched_barcodes = [None] * len(read_seqs)

            with ProcessPoolExecutor(max_workers=threads, initializer=init_barcode_match_worker, initargs=(barcode_seqs, barcode_seqs_full)) as executor:
                for start, batch_min, batch_matches in executor.map(match_barcode_batch, batches):
                    min_dist[start:start + len(batch_min)] = batch_min
                    for i, m in enumerate(batch_matches):
                        matched_barcodes[start + i] = m

            result = pd.DataFrame({
                'sequence': read_seqs,
                'min_hamming_distance': min_dist,
                'barcode_matches': matched_barcodes,
            })
            result.to_csv(output.matches, index=False)


rule apply_barcode_matches_quality_table_to_tile_zscored:
    """ Joins the well-wide barcode match lookup back onto one tile's corrected read table.
    This is a cheap dict-based .map(), not a repeat alignment -- see the discussion in
    cell_table_attempts.ipynb on why this stays fast even at full-well scale.
    """
    input:
        table = sequencing_dir + '{well}_grid{grid_size}/{tile}/{segmentation_type}_quality{params}.csv',
        matches =  sequencing_dir + '{well}_grid{grid_size}/{segmentation_type}_zscore_barcode_matches{params}.csv',
    output:
        table = sequencing_dir + '{well}_grid{grid_size}/{tile}/{segmentation_type}_zscored_matched_quality{params}.csv',
    wildcard_constraints:
        grid_size = '\d+',
        params = params_regex('min', 'max', 'num', ),
        tile = 'tile\d+x\d+y',
    resources:
        mem_mb = lambda wildcards, input: 5000 + size_mb(input) * 3
    run:
        import pandas as pd
        import starcall.reads

        table = pd.read_csv(input.table, index_col=0)
        table['max_seq'] = table.reads.sequences
        matches = pd.read_csv(input.matches)
        dist_lookup = dict(zip(matches['sequence'], matches['min_hamming_distance']))
        match_lookup = dict(zip(matches['sequence'], matches['barcode_matches']))

        table['min_hamming_distance'] = table['max_seq'].map(dist_lookup)
        table['barcode_matches'] = table['max_seq'].map(match_lookup)

        table.to_csv(output.table)


rule apply_barcode_matches_quality_table_to_fullwell_zscored:
    """ Joins the well-wide barcode match lookup back onto one tile's corrected read table.
    This is a cheap dict-based .map(), not a repeat alignment -- see the discussion in
    cell_table_attempts.ipynb on why this stays fast even at full-well scale.
    """
    input:
        table = sequencing_dir + '{well}/{segmentation_type}_quality{params}.csv',
        matches =  sequencing_dir + '{well}/{segmentation_type}_zscore_barcode_matches{params}.csv',
    output:
        table = sequencing_dir + '{well}/{segmentation_type}_zscored_matched_quality{params}.csv',
    wildcard_constraints:
        params = params_regex('min', 'max', 'num', ),
    resources:
        mem_mb = lambda wildcards, input: 5000 + size_mb(input) * 3
    run:
        import pandas as pd
        import starcall.reads

        table = pd.read_csv(input.table, index_col=0)
        table['max_seq'] = table.reads.sequences
        matches = pd.read_csv(input.matches)
        dist_lookup = dict(zip(matches['sequence'], matches['min_hamming_distance']))
        match_lookup = dict(zip(matches['sequence'], matches['barcode_matches']))

        table['min_hamming_distance'] = table['max_seq'].map(dist_lookup)
        table['barcode_matches'] = table['max_seq'].map(match_lookup)

        table.to_csv(output.table)

rule combine_reads_in_cells_with_scores_zscored_quality_version:
    """
    All reads in each cell are sorted by the read count.
    """
    input:
        #sequencing_dir + '{well}_grid{grid_size}/{tile}/{segmentation_type}{approach}_corrected_matched_{method}{params}.csv',
        table = sequencing_dir + '{path}/{segmentation_type}_zscored_matched_quality{params}.csv',
        cell_table = segmentation_dir + '{path}/{segmentation_type}.csv',
    output:
        table = sequencing_dir + '{path}/{segmentation_type}_reads{params}.csv',
    wildcard_constraints:
        maxreads = '|_maxreads\d+',
        params = params_regex('min', 'max', 'num', ),
    resources:
        mem_mb = lambda wildcards, input: 5000 +  size_mb(input) * 50
    run:
        import pandas
        import numpy as np
        import starcall.utils
        import starcall.reads

        table = pandas.read_csv(input.table, index_col=0)
        table = table.loc[table['cell']!=0,:]
        
        #prep to match empty cell rows
        cell_table = pandas.read_csv(input.cell_table, index_col=0)
        full_cells_index = range(1, len(cell_table.index) + 1)

        #group into cells
        #sort by phred score, so highest phred score is first for reads with multiple dots
        cell_col='cell'
        seq_col='max_seq'
        phred_col='min_chastity' #change to chastity values 
        barcode_col = 'barcode_matches'
        hamming_col = 'min_hamming_distance'
        sub = table[[cell_col, seq_col, barcode_col, hamming_col, phred_col]]

        sub = sub.sort_values(phred_col, ascending = False)
        grouped = (
            sub.groupby([cell_col, seq_col, barcode_col, hamming_col])[phred_col]
            .apply(lambda x: ':'.join(x.round(4).astype(str)))
            .reset_index(name='chastities')
        ) #keep only 4 decimal places for the final table
        grouped['count'] = sub.groupby([cell_col, seq_col]).size().values
        grouped['first_val'] = grouped.chastities.apply(lambda x: float(x.split(':')[0]))

        # order each cell's unique seqs by how many reads support them (most first)
        grouped = grouped.sort_values([cell_col, 'count', 'first_val'], ascending=[True, False, False])
        grouped['seq_rank'] = grouped.groupby(cell_col).cumcount() 
        debug (grouped.seq_rank.max())

        seq_wide = grouped.pivot(index=cell_col, columns='seq_rank', values=seq_col)
        seq_wide.columns = [f'read_{i}' for i in seq_wide.columns]
        #count unique reads per cell
        seq_wide['num_reads'] = seq_wide.count(axis=1) #count non-nan entries per row

        quality_wide = grouped.pivot(index=cell_col, columns='seq_rank', values='chastities')
        quality_wide.columns = [f'chastities_{i}' for i in quality_wide.columns]

        barcodes_wide = grouped.pivot(index=cell_col, columns='seq_rank', values='barcode_matches')
        barcodes_wide.columns = [f'barcode_matches_{i}' for i in barcodes_wide.columns]

        hdist_wide = grouped.pivot(index=cell_col, columns='seq_rank', values='min_hamming_distance')
        hdist_wide.columns = [f'barcode_hamming_dist_{i}' for i in hdist_wide.columns]

        count_wide = grouped.pivot(index=cell_col, columns='seq_rank', values='count')
        count_wide.columns = [f'count_{i}' for i in count_wide.columns]
        #setup total_reads col for cell (sum of occurences of each unique read)
        count_wide = count_wide.fillna(0)
        count_wide['total_count'] = count_wide[count_wide.columns].sum(axis=1)

        ordered_cols = [c for i in range(0, seq_wide.shape[1] -1 ) for c in (f'read_{i}', f'count_{i}', f'chastities_{i}', f'barcode_matches_{i}', f'barcode_hamming_dist_{i}')]
        ordered_cols = ['num_reads', 'total_count', 'dot_indicies'] + ordered_cols

        #':'-joined indices of every dot in the cell, pointing back to rows of this tile's quality table (used by the viewer)
        dot_indicies = table.index.to_series().groupby(table[cell_col]).agg(lambda x: ':'.join(map(str, x))).rename('dot_indicies')

        table = pd.concat([seq_wide, count_wide, dot_indicies, quality_wide, barcodes_wide, hdist_wide], axis=1)[ordered_cols]

        #insert empty rows for any cells with no reads 
        table = table.reindex(full_cells_index)
        cell_reads = table.set_index(cell_table.index)
        cell_reads.to_csv(output.table)


#barcode matching approaches - 3 main approaches, need to choose a final winner and give it the expected barcode values? 
#what to merge in from the barcodes table....
#check merge_final_table in sequencing....




def mult_zero_match_filter(row):
    #return true if more than one barcode_hamming_dist_{i} column is equal to 0 
    zero_count = sum(1 for ax in row.axes[0] if ax.startswith('barcode_hamming_dist_') and row[ax] == 0)
    return zero_count > 1

def exact_zero_match_filter(row):
    #return true if more than one barcode_hamming_dist_{i} column is equal to 0 
    zero_count = sum(1 for ax in row.axes[0] if ax.startswith('barcode_hamming_dist_') and row[ax] == 0)
    return zero_count == 1


def get_matched_corrected_tables_zscored(wildcards):
    grid_size = int(wildcards.grid_size)
    numbers = ['{:02}'.format(i) for i in range(grid_size)]
    #'{path}/{segmentation_type}_zscored_matched_cell_quality{params}.csv'
    return expand(sequencing_dir + '{well}_grid{grid_size}/tile{x}x{y}y/{segmentation_type}_zscored_matched_cell_quality{params}.csv', x=numbers, y=numbers, allow_missing=True)

rule calculate_per_tile_stats_zscored: 
    input:
        corrected_tables = get_matched_corrected_tables_zscored,
    output:
        summary_csv = sequencing_dir + '{well}_grid{grid_size}/{segmentation_type}_zscored_per_cell_summaries{params}.csv',
    wildcard_constraints:
        params = params_regex('min', 'max', 'num', ),
    resources:
        mem_mb = lambda wildcards, input: 5000 + size_mb(input) * 2
    run:
        import pandas as pd

        stats_df = pd.DataFrame({'tile':[], 'frac_single_reads_only':[], 'frac_one_exact_match':[], 'frac_mult_exact_matches':[], 'frac_top_read_exact_match':[]})
        summary_rows = []
        for table_path in input.corrected_tables:
            tile_match = re.search(r'/(tile\d+x\d+y)/', table_path)
            tile = tile_match.group(1) if tile_match else table_path

            grouped_table = pd.read_csv(table_path, index_col = 0)

            grouped_table['mult_zero_matches'] = grouped_table.apply(lambda row: mult_zero_match_filter(row), axis = 1)
            grouped_table['exact_zero_matches'] = grouped_table.apply(lambda row: exact_zero_match_filter(row), axis = 1)

            frac_single_reads_only = grouped_table[(grouped_table.num_reads == grouped_table.total_count)].shape[0]/grouped_table.shape[0]
            frac_has_zero_match = grouped_table[grouped_table.exact_zero_matches].shape[0]/grouped_table.shape[0]
            frac_has_mult_zero_matches = grouped_table[grouped_table.mult_zero_matches].shape[0]/grouped_table.shape[0]
            frac_top_read_edit_dist_0 = grouped_table[(grouped_table.barcode_hamming_dist_0 == 0) & (grouped_table.count_0 != 1)].shape[0]/grouped_table.shape[0]

            summary_rows.append({
                'tile':tile,
                'frac_single_reads_only':frac_single_reads_only,
                'frac_has_one_zero_match':frac_has_zero_match,
                'frac_has_mult_zero_matches':frac_has_mult_zero_matches,
                'frac_top_reads_edit_dist_0':frac_top_read_edit_dist_0,
            })
        summary_table = pd.DataFrame(summary_rows)
        summary_table.to_csv(output.summary_csv)
