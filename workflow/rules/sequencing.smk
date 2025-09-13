import os
import glob
import re


rule find_dots:
    input:
        stitching_dir + '{path}/raw.tif'
    output:
        sequencing_dir + '{path}/bases.csv',
        #sequencing_dir + '{path}/dot_filter.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 10 + 15000
    threads: 4
    run:
        import numpy as np
        import tifffile
        import starcall.dotdetection
        import starcall.correction
        import skimage.morphology

        min_sigma = 1#2
        max_sigma = 2#3
        num_sigma = 7

        full_well = tifffile.memmap(input[0], mode='r')
        image = full_well[...,2:,:,:].astype(np.float32, copy=True)
        del full_well

        if np.all(image == 0):
            reads = ReadSet()
        else:
            #dot_filter = starcall.dotdetection.dot_filter_new(image)
            #tifffile.imwrite(output[1], dot_filter)
            reads = starcall.dotdetection.detect_dots(
                image,
                min_sigma = min_sigma,
                max_sigma = max_sigma,
                num_sigma = num_sigma,
                copy = False,
            )

        reads.to_csv(output[0])


rule call_raw_reads:
    input:
        bases = sequencing_dir + '{path}/bases.csv',
        cells = segmentation_dir + '{path}/{segmentation_type}_mask_downscaled.tif',
    output:
        table = sequencing_dir + '{path}/{segmentation_type}_raw_reads.csv',
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 10
    run:
        import tifffile
        import numpy as np
        import pandas
        import csv
        import matplotlib.pyplot as plt
        import starcall.reads

        cells = tifffile.imread(input.cells)
        table = pandas.read_csv(input.bases, index_col=0)

        xposes, yposes = np.round(table.reads.position.T).astype(int)
        table['cell'] = cells[xposes,yposes]

        table.to_csv(output.table)


rule calculate_distance_matrix:
    input:
        bases = sequencing_dir + '{path}/bases.csv',
        cells = segmentation_dir + '{path}/{segmentation_type}_mask_downscaled.tif',
    output:
        table = sequencing_dir + '{path}/{segmentation_type}_reads_distance_matrix.csv'
    resources:
        mem_mb = lambda wildcards, input: 15000 + input.size_mb * 10
    run:
        import tifffile
        import numpy as np
        import pandas
        import starcall.reads
        import csv

        normalization = 'none'
        positional_weight = 100
        value_weight = 0
        sequence_weight = 1

        table = pandas.read_csv(input.bases, index_col=0)
        cells = tifffile.imread(input.cells)

        if normalization != 'none':
            table.reads.normalize(method=normalization)

        distance_matrix = starcall.reads.distance_matrix(
            table, cells=cells,
            distance_cutoff=50,
            positional_weight=positional_weight,
            value_weight=value_weight,
            sequence_weight=sequence_weight,
            debug=True, progress=True,
        )

        with open(output.table, 'w') as ofile:
            writer = csv.DictWriter(ofile, ['i', 'j', 'distance'])
            writer.writeheader()

            for pair, dist in distance_matrix.items():
                writer.writerow(dict(i=pair[0], j=pair[1], distance=dist))


rule cluster_reads:
    input:
        distances = sequencing_dir + '{path}/{segmentation_type}_reads_distance_matrix.csv'
    output:
        clusters = sequencing_dir + '{path}/{segmentation_type}_reads_clusters.csv',
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 10
    run:
        import numpy as np
        import csv
        import starcall.reads

        threshold = 0.5
        linkage = 'min'

        distance_matrix = {}
        with open(input.distances) as ifile:
            reader = csv.DictReader(ifile)
            for row in reader:
                i, j, distance = int(row['i']), int(row['j']), float(row['distance'])
                distance_matrix[i,j] = distance

        cluster_indices = starcall.reads.cluster_reads(
            distance_matrix,
            threshold=threshold,
            linkage=linkage,
            debug=True, progress=True,
        )

        with open(output.clusters, 'w') as ofile:
            ofile.write('cluster\n')
            ofile.write(''.join(str(cluster) + '\n' for cluster in cluster_indices))



rule combine_reads:
    input:
        raw_reads = sequencing_dir + '{path}/{segmentation_type}_raw_reads.csv',
        clusters = sequencing_dir + '{path}/{segmentation_type}_reads_clusters.csv',
    output:
        table = sequencing_dir + '{path}/{segmentation_type}_clustered_reads.csv',
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 25
    run:
        import tifffile
        import numpy as np
        import pandas
        import csv
        import matplotlib.pyplot as plt
        import starcall.reads

        table = pandas.read_csv(input.raw_reads, index_col=0)
        clusters = np.loadtxt(input.clusters, skiprows=1, dtype=int, delimiter=',').reshape(-1)
        table['cluster'] = clusters
        table['count'] = np.ones(len(clusters))

        combined = table.groupby('cluster')
        combined = combined.agg(**table.reads.aggfuncs(position='mean', values='sum', count='sum', cell=pandas.Series.mode))
        combined.reads.normalize()

        combined.to_csv(output.table)


rule combine_cell_reads:
    input:
        table = sequencing_dir + '{path}/{segmentation_type}_clustered_reads.csv',
    output:
        table = sequencing_dir + '{path}/{segmentation_type}_reads_partial.csv',
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 50
    run:
        import pandas
        import numpy as np
        import starcall.utils
        import starcall.reads

        max_reads = 2

        table = pandas.read_csv(input.table)
        table = table.loc[table['cell']!=0,:]

        cell_reads = table.sort_values(['cell', 'count'], ascending=False).groupby('cell')
        cell_reads = cell_reads.head(max_reads)
        cell_reads = cell_reads.reads.to_cell_table()
        #read_table = cell_reads.head(max_reads).to_table(columns=['cell', 'read', 'count', 'quality', 'read_index'], sequences=True, qualities=True)
        #read_table = read_table.set_index('cell')
        #read_table['total_count'] = [read_set.attrs['count'].sum() for read_set in cell_reads]
        cell_reads.to_csv(output.table)


#ruleorder: merge_grid > find_dots
#ruleorder: merge_grid > segment_cells
#ruleorder: link_grid > find_dots
#ruleorder: link_grid > segment_cells
#ruleorder: merge_grid > segment_cells_bases
ruleorder: segment_cells > segment_cells_bases


rule annotate_dots:
    input:
        image = stitching_dir + '{path}/raw.tif',
        bases = sequencing_dir + '{path}/bases.csv',
        #'tmp_dot_greyimage.tif',
        #sequencing_dir + '{path}/clusters.csv',
    output:
        qc_dir + '{path}/annotated.tif',
        #qc_dir + '{path}/annotated_clusters.tif',
        #qc_dir + '{path}/annotated_grey.tif',
    resources:
        mem_mb = 25000
    run:
        import numpy as np
        import starcall.utils
        import tifffile
        import skimage.draw
        import starcall.reads
        import pandas

        reads = starcall.reads.ReadSet.from_table(pandas.read_csv(input.bases, index_col=0))
        #image = tifffile.imread(input[0])
        image = tifffile.memmap(input[0], mode='r')[0]
        #poses = np.loadtxt(input[1], delimiter=',')[:,:2].astype(int)
        #clusters = open(input[3]).read().strip().split()[1:]
        #clusters = np.array([int(cluster) for cluster in clusters])

        marked_image = starcall.utils.mark_dots(image, reads.positions.astype(int))
        tifffile.imwrite(output[0], marked_image)

        """
        marked_image[-1] = 0
        for cluster in range(clusters.max()+1):
            cluster_poses = poses[clusters==cluster]
            center = cluster_poses.mean(axis=0).astype(int)
            for pos in cluster_poses:
                line = skimage.draw.line(*center, *pos)
                marked_image[-1,line[0],line[1]] = 1

        tifffile.imwrite(output[1], marked_image)

        greyimage = tifffile.imread(input[2])
        tifffile.imwrite(output[2], starcall.utils.mark_dots(greyimage[None,:,:], poses))
        """


#ruleorder: tabulate_cells > merge_grid

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

rule merge_final_tables:
    input:
        cell_table = sequencing_dir + '{path}/{segmentation_type}_reads_partial.csv',
        aux_data = get_aux_data,
    output:
        full_table = sequencing_dir + '{path}/{segmentation_type}_reads.csv',
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 250
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
            barcodes = []
            num_cycles = len(cell_table['read_0'].iloc[0])
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

            #library = starcall.sequencing.BarcodeLibrary(barcodes)

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

            new_rows = [pandas.Series(dtype=object) if i == -1 else aux_table.iloc[i,:] for i in indices[:,0]]
            #new_rows = [pandas.Series(dtype=object) if (i == -1 and not multiple) else aux_table.iloc[i,:] for i, multiple in zip(indices[:,0], multiple_matches)]
            debug ('making new table', len(new_rows))
            new_table = pandas.DataFrame(new_rows, index=cell_table.index)
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

#ruleorder: call_reads > merge_grid
#ruleorder: match_masks > merge_grid
#ruleorder: merge_final_tables > merge_grid
#ruleorder: calc_features > merge_grid


##################################################
## Merging read tables
##################################################

def get_grid_filenames_seq(wildcards):
    grid_size = int(wildcards.grid_size)
    numbers = ['{:02}'.format(i) for i in range(grid_size)]
    return expand(sequencing_dir + '{well}_grid{grid_size}/tile{x}x{y}y/{segmentation_type}_reads.csv', x=numbers, y=numbers, allow_missing=True)

rule merge_grid_read_tables:
    input:
        tables = get_grid_filenames_seq,
        #composite = stitching_dir + '{well}_grid{grid_size}/grid_composite.json',
    output:
        table = sequencing_dir + '{well}_grid{grid_size,\d+}/{segmentation_type}_reads.csv',
    resources:
        #mem_mb = lambda wildcards, input: input.size_mb * 50 + 5000
        mem_mb = 5000
    run:
        import constitch

        #composite = constitch.load(input.composite)

        def row_func(row):
            i = row['file_index']
            #box = composite.boxes[i]
            row['tile_index'] = i
            row['tile_x'] = i // int(wildcards.grid_size)
            row['tile_y'] = i % int(wildcards.grid_size)

        merge_csv_files(input.tables, output.table, extra_columns=['tile_x', 'tile_y', 'tile_index'], row_func=row_func)


if config.get('sequencing_grid_size', 1) != 1:
    rule link_merged_grid_reads:
        input:
            sequencing_dir + '{path_nogrid}_grid' + str(config['sequencing_grid_size']) + '/{segmentation_type}_reads.csv',
        output:
            sequencing_dir + '{path_nogrid}/{segmentation_type}_reads.csv',
        localrule: True
        wildcard_constraints:
            path_nogrid = '((?!_grid)[^.])*',
        shell:
            "cp -l {input[0]} {output[0]}"

    ruleorder: link_merged_grid_reads > merge_final_tables

