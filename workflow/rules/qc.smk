import os
import sys
import glob
import time

##################################################
## Quality control plots
##################################################

rule make_qc_read_plots:
    input:
        full_table = sequencing_dir + '{path}/{segmentation_type}_reads.csv',
        barcodes = get_aux_data,
    output:
        plot = qc_dir + '{path}/{segmentation_type}_reads.svg',
        plots = [qc_dir + '{path}/{segmentation_type}_reads_plot' + str(i) + '.svg' for i in range(16)],
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 10
    run:
        import pandas
        import numpy as np
        import starcall.utils
        import matplotlib
        import matplotlib.pyplot as plt

        debug (matplotlib.rcParams['font.size'])
        matplotlib.rcParams.update({'xtick.labelsize': 10, 'ytick.labelsize': 10, 'font.size': 12})
        #12 for ticks, 14 for labels, important to be consistent

        read_table = pandas.read_csv(input.full_table, index_col=0)
        num_cycles = len(read_table['read_0'].iloc[0])
        library_paths = [path for path in input.barcodes if path.count('barcodes.csv')]
        has_library = len(library_paths) != 0
        if has_library:
            library_table = pandas.read_csv(library_paths[0])
            barcodes = [(barc + '???????????????????')[:num_cycles] for barc in library_table[library_table.columns[0]]]
            barcode_mat = np.array(list(map(list, barcodes)))
            debug (barcode_mat.shape, 'barcode mat')

        if 'num_reads' in read_table.columns:
            debug ('normal loading')
            reads = starcall.utils.read_multicolumn(read_table, 'read', length_column='num_reads')
            counts = starcall.utils.read_multicolumn(read_table, 'count', length_column='num_reads')
            #qualities = starcall.utils.read_multicolumn(read_table, 'quality', length_column='num_reads')
        else:
            debug ('other loading')
            reads_mat = starcall.utils.read_multicolumn(read_table, 'read')
            counts_mat = starcall.utils.read_multicolumn(read_table, 'count')
            #qualities_mat = starcall.utils.read_multicolumn(read_table, 'quality')
            reads, counts, qualities = [], [], []

            for i in range(len(reads_mat)):
                reads.append([read for read, count in zip(reads_mat[i], counts_mat[i]) if count != 0])
                counts.append([count for count in counts_mat[i] if count != 0])
                #qualities.append([quality for quality, count in zip(qualities_mat[i], counts_mat[i]) if count != 0])

        debug (counts[:10])
        for i in range(len(reads)):
            num_reads = min(len(counts[i]), len(reads[i]))
            #num_reads = min(len(counts[i]), len(reads[i]), len(qualities[i]))
            reads[i] = [reads[i][j] for j in range(num_reads) if counts[i][j] != 0]
            #qualities[i] = [qualities[i][j] for j in range(num_reads) if counts[i][j] != 0]
            tmp = [counts[i][j] for j in range(num_reads) if counts[i][j] != 0]
            counts[i] = tmp
        debug (counts[:10])

        if 'total_count' not in read_table.columns:
            #read_table['total_count'] = read_table['num_reads']
            read_table['total_count'] = list(map(sum, counts))

        double_barcode = 'matched_read_index_1' in read_table.columns



        def make_plots(fig, axes):
            read_mat = []
            for read_set in reads:
                read_mat.extend(map(list, read_set))

            debug (len(reads), len(read_mat))
            read_mat = np.array(read_mat)
            debug (read_mat.shape)

            for i,let in enumerate('GTAC'):
                if has_library:
                    axes[0,0].plot(np.mean(barcode_mat == let, axis=0), '--C' + str(i))
                    debug ('axis.plot({}, "--C{}", label="{}")'.format(np.mean(barcode_mat == let, axis=0).tolist(), i, let))
                axes[0,0].plot(np.mean(read_mat == let, axis=0), '-C' + str(i), label=let)
                debug ('axis.plot({}, "-C{}", label="{}")'.format(np.mean(read_mat == let, axis=0).tolist(), i, let))

            #axes[0,0].set_title('Nucleotide frequencies' + ('\n(dashed is library)' if has_library else ''))
            axes[0,0].set_ylabel('Nucleotide frequencies' + ('\n(dashed is library)' if has_library else ''))
            axes[0,0].set_xlabel('Cycle')
            axes[0,0].set_ylim(ymin=0)
            axes[0,0].legend(loc='upper left')
            #if 'blainey' in wildcards.path:
                #axes[0,0].set_title('Image set 1.1: Feldman method')
            #else:
                #axes[0,0].set_title('Image set 1.1: STARCall')

            total_cells = max(read_table.index)
            labels, values = np.unique(read_table['total_count'], return_counts=True)
            labels, values = labels[labels<10], values[labels<10]
            labels = np.array([0] + list(labels))
            values = np.array([total_cells - len(read_table.index)] + list(values))
            debug ('labels = ', labels[:30].tolist())
            debug ('values = ', values[:30].tolist())
            axes[0,1].bar(labels[:30], values[:30], width=1)
            #axes[0,1].set_title('Read count of cells')
            axes[0,1].set_ylabel('Count')
            axes[0,1].set_xlabel('Total read count per cell' if double_barcode else 'Read count per cell') # total read count for double barcode
            axes[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            if 'blainey' in wildcards.path:
                axes[0,1].set_title('Image set 2.1: Feldman method')
            else:
                axes[0,1].set_title('Image set 2.1: STARCall')
            # todo make heatmap of first and second barcode edit distance in double barcode
            # total read count, combined edit distance
            # edit distance to first matching barcode pair
            # expected combination in every cell
            # permute in cells, show 

            # Final plots:
            # lmna: total read count, edit distance hist
            # pten: stacked total read count for ipsc and neuron, stacked edit distance, cell level heatmap

            axes[0,2].bar(*np.unique(list(map(len, reads)), return_counts=True), width=1)
            #axes[0,2].bar(*np.unique(read_table['num_reads'], return_counts=True), width=1)
            #axes[0,2].set_title('Number of unique reads per cell')
            axes[0,2].set_ylabel('Count')
            axes[0,2].set_xlabel('Unique read count')

            if has_library:
                labels, values = np.unique(read_table['edit_distance'], return_counts=True)
                #debug (labels)
                #none_vals = values[0]
                total_count = values.sum()
                #labels = ['None'] + list(map(str, labels))
                #values = [0] + list(values)
                labels = list(map(str, labels))
                values = list(values)
                #debug (labels)
                axes[0,3].bar(labels, values, width=1, color='grey', label='mul.')
                labels, values = np.unique(read_table.loc[~pandas.isna(read_table['aaChanges']),'edit_distance'], return_counts=True)
                assert -1 not in labels

                total_not_included = total_count - values.sum()
                labels = list(map(str, labels))
                values = list(values)
                #labels = ['None'] + list(map(str, labels))
                #values = [total_not_included] + list(values)
                axes[0,3].bar(labels, values, width=1, label='sing.')

                #plotting all reads that were removed in the none bar
                #axes[0,3].bar(['None'], [total_not_included], width=1)

                #axes[0,3].set_title('Edit distance, single v multiple')
                axes[0,3].legend()
                axes[0,3].set_xlabel('Combined edit distance' if double_barcode else 'Edit distance')
                axes[0,3].set_ylabel('Count')

                matched_reads = []
                matched_barcs = []
                for i in read_table.index:
                    matched_index = read_table.loc[i,'matched_read_index_0']
                    if matched_index != -1:
                        matched_reads.append(read_table.loc[i,'read_{}'.format(matched_index)])
                        matched_barcs.append((read_table.loc[i,'matched_barcode_0'] + '????????????')[:num_cycles])

                matched_reads = np.array(list(map(list, matched_reads)))
                matched_barcs = np.array(list(map(list, matched_barcs)))

                debug (matched_reads.shape, 'matched_reads')
                debug (matched_barcs.shape, 'matched_barcs')
                debug ((matched_reads != matched_barcs).shape, 'matched_barcs')

                tmp_edit_dist = np.sum(matched_reads != matched_barcs, axis=1)
                print (tmp_edit_dist)
                print (np.unique(tmp_edit_dist, return_counts=True))
                #ksdjf

                axes[1,0].set_title('Errors between reads and library over cycles')
                axes[1,0].plot(np.sum(matched_reads != matched_barcs, axis=0))

                axes[1,1].set_title('Errors by read nucleotide')
                axes[1,2].set_title('Errors by library nucleotide')
                axes[1,3].set_title('Errors by library nucleotide to read nucleotide')

                for let in 'GTAC':
                    axes[1,1].plot(np.sum((matched_reads != matched_barcs) & (matched_reads == let), axis=0), label=let)
                    axes[1,2].plot(np.sum((matched_reads != matched_barcs) & (matched_barcs == let), axis=0), label=let)
                    for let2 in 'GTAC':
                        axes[1,3].plot(np.sum((matched_reads != matched_barcs)
                                & (matched_barcs == let) & (matched_reads == let2), axis=0),
                                label=let + '->' + let2)

                axes[1,0].set_ylim(ymin=0)
                axes[1,1].set_ylim(ymin=0)
                axes[1,2].set_ylim(ymin=0)
                axes[1,3].set_ylim(ymin=0)

                axes[1,1].legend()
                axes[1,2].legend()
                axes[1,3].legend()

                #coverage, e_values = starcall.sequencing.calculate_library_coverage
                #axes[2,0].plot(

                max_dist = max(read_table['edit_distance'].max(), read_table['edit_distance_2'].max())
                max_dist = 5 if double_barcode else 4
                values = np.zeros((max_dist + 1, max_dist + 1))
                for dist1, dist2 in zip(read_table['edit_distance'], read_table['edit_distance_2']):
                    if dist1 < 0 or dist2 < 0: continue
                    values[min(dist1, max_dist),min(dist2, max_dist)] += 1

                total = int(values.sum())
                max_val = values.max()
                values[values==0] = np.nan

                imshow = axes[2,0].imshow(values)
                #fig.colorbar(imshow, ax=axes[2,0])

                for x in range(values.shape[0]):
                    for y in range(values.shape[1]):
                        if np.isnan(values[x,y]): continue
                        color = 'white' if values[x,y] / max_val < 0.66 else 'black'
                        axes[2,0].text(y, x, '{:.1f}%'.format(values[x,y] / total * 100), ha='center', va='center', color=color, size=10)

                axes[2,0].set_xticks(list(range(max_dist + 1)))
                axes[2,0].set_yticks(list(range(max_dist + 1)))
                axes[2,0].set_xticklabels(list(map(str, range(max_dist))) + [str(max_dist) + '+'])
                axes[2,0].set_yticklabels(list(map(str, range(max_dist))) + [str(max_dist) + '+'])
                if not double_barcode:
                    axes[2,0].set_ylabel('Edit distance to best\nmatching barcode')
                    axes[2,0].set_xlabel('Edit distance to second best\nmatching barcode')
                else:
                    axes[2,0].set_ylabel('Total edit distance to best\nmatching barcode pair')
                    axes[2,0].set_xlabel('Total edit distance to second best\nmatching barcode pair')
                debug (axes[2,0].get_xticklabels())
                #axes[2,0].set_title('Edit distance to first and second matching barcode')

                """
                matched_qualities = []
                unmatched_qualities = []
                for i in range(len(qualities)):
                    if read_table['edit_distance'].iloc[i] == 0:
                        index = read_table['matched_read_index_0'].iloc[i]
                        for j in range(len(qualities[i])):
                            if j == index:
                                for k in range(int(counts[i][j])):
                                    matched_qualities.append(qualities[i][j])
                            else:
                                for k in range(int(counts[i][j])):
                                    unmatched_qualities.append(qualities[i][j])

                        #matched_qualities.append(qualities[i,index])
                        #unmatched_qualities.extend(qualities[i,:index])
                        #unmatched_qualities.extend(qualities[i,index+1:])

                debug (matched_qualities[:50])
                debug (unmatched_qualities[:50])

                matched_qualities = np.sort(matched_qualities)
                unmatched_qualities = np.sort(unmatched_qualities)
                debug (matched_qualities.shape, unmatched_qualities.shape)

                points = np.linspace(unmatched_qualities.min(), matched_qualities.max(), 15)
                above_matched = np.array([matched_qualities.shape[0] - np.searchsorted(matched_qualities, point) for point in points])
                above_unmatched = np.array([unmatched_qualities.shape[0] - np.searchsorted(unmatched_qualities, point) for point in points])

                axes[2,1].plot(points, (above_matched + above_unmatched) / (len(matched_qualities) + len(unmatched_qualities)), label='Above thresh')
                axes[2,1].plot(points, above_matched / (above_matched + above_unmatched), label='Barcode match')

                axes[2,1].legend()
                axes[2,1].set_xlabel('Read quality threshold')
                axes[2,1].set_ylabel('Reads (%)')
                """

                matched_read_counts = []
                for i in range(len(reads)):
                    index = read_table['matched_read_index_0'].iloc[i]
                    matched_read_counts.append(counts[i][index] if index != -1 else 0)

                labels, values = np.unique(matched_read_counts, return_counts=True)
                if len(values) > 5:
                    values[5] = values[5:].sum()
                    labels = list(map(str, labels[:5].astype(int))) + ['5+']
                    values = values[:6]

                axes[2,2].bar(labels, values, width=1)
                #axes[2,2].set_title('Read count of cells (exact barcode match)')
                axes[2,2].set_ylabel('Count')
                axes[2,2].set_xlabel('Barcode count of cells')


                values = np.zeros((5, 5))

                for i in range(len(reads)):
                    index = read_table['matched_read_index_0'].iloc[i]
                    if index != -1:
                        total = int(sum(counts[i]))
                        values[min(int(counts[i][index]) - 1, 4), min(total - 1, 4)] += 1

                total = int(values.sum())
                max_val = values.max()
                values[values==0] = np.nan

                imshow = axes[2,3].imshow(values)

                for x in range(values.shape[0]):
                    for y in range(values.shape[1]):
                        if np.isnan(values[x,y]): continue
                        color = 'white' if values[x,y] / max_val < 0.66 else 'black'
                        axes[2,3].text(y, x, '{:.2f}%'.format(values[x,y] / total * 100), ha='center', va='center', color=color)

                axes[2,3].set_xlabel('Total reads')
                axes[2,3].set_ylabel('Top barcode reads')

                axes[2,3].set_xticks(list(range(5)))
                axes[2,3].set_yticks(list(range(5)))
                axes[2,3].set_xticklabels(list(map(str, range(1, 5))) + ['5+'])
                axes[2,3].set_yticklabels(list(map(str, range(1, 5))) + ['5+'])


                if 'matched_read_index_1' in read_table.columns:
                    max_read_count = 6
                    values = np.zeros((max_read_count + 1, max_read_count + 1))

                    for i in range(len(reads)):
                        #count1 = int(counts[i][0])
                        #count2 = 0 if len(counts[i]) <= 1 else int(counts[i][1])
                        #values[:min(count1,max_read_count)+1,:min(count2,max_read_count)+1] += 1
                        index1, index2 = read_table['matched_read_index_0'].iloc[i], read_table['matched_read_index_1'].iloc[i]
                        if index1 != -1 and index2 != -1:
                            count1, count2 = int(counts[i][index1]), int(counts[i][index2])
                            if count2 < count1:
                                count1, count2 = count2, count1
                            #values[:min(count1,max_read_count)+1,:min(count2,max_read_count)+1] += 1
                            values[min(count1,max_read_count),min(count2,max_read_count)] += 1

                    #max_val = values[0,0]
                    #total = values[0,0]
                    max_val = int(values.max())
                    total = int(values.sum())

                    for i in range(max_read_count + 1):
                        for j in range(i):
                            values[i,j] = np.nan

                    imshow = axes[3,0].imshow(values)

                    for x in range(values.shape[0]):
                        for y in range(values.shape[1]):
                            debug (values[x,y], 'value')
                            if np.isnan(values[x,y]): continue
                            color = 'white' if values[x,y] / max_val < 0.66 else 'black'
                            axes[3,0].text(y, x, '{:.1f}%'.format(values[x,y] / total * 100), ha='center', va='center', color=color)


                total_count = 0
                errors = np.zeros((num_cycles, 4, 4), dtype=int)
                for i in range(len(reads)):
                    matched_read_index = read_table['matched_read_index_0'].iloc[i]
                    if read_table['edit_distance'].iloc[i] <= 1 and matched_read_index != -1:
                        barcode = read_table['matched_barcode_0'].iloc[i]
                        read = reads[i][matched_read_index]
                        count = counts[i][matched_read_index]
                        #if count < 6: continue
                        diffs = [let1 != let2 for let1, let2 in zip(read, barcode)]
                        if sum(diffs) != 1: continue
                        index = diffs.index(True)
                        read_index, barc_index = 'GTAC'.index(read[index]), 'GTAC'.index(barcode[index])
                        errors[index,read_index,barc_index] += count
                        total_count += count
                        #dists = [dist for dist in dists if dist == 1]
                        #edit_distance = sum(dist * count for dist, count in zip(dists, counts[i]))
                        #total_count += num_cycles * sum(counts[i])

                #axes[3,0].set_title('error rate {}'.format(errors / total_count * num_cycles))
                axes[3,0].plot(errors.sum(axis=(1,2)) / total_count)
                axes[3,0].set_xlabel('Cycle')
                axes[3,0].set_ylabel('Error rate')


                matched_edit_distances = []
                total_count = 0
                #errors = [[] for i in range(num_cycles)]
                errors = np.zeros((num_cycles, 4, 4), dtype=int)
                for i in range(len(reads)):
                    if read_table['edit_distance'].iloc[i] == 0 and read_table['matched_read_index_0'].iloc[i] != -1:
                        barcode = read_table['matched_barcode_0'].iloc[i]
                        dists = [sum(let1 != let2 for let1, let2 in zip(read, barcode)) for read in reads[i]]
                        for read, count in zip(reads[i], counts[i]):
                            diffs = [let1 != let2 for let1, let2 in zip(read, barcode)]
                            if sum(diffs) != 1: continue
                            index = diffs.index(True)
                            read_index, barc_index = 'GTAC'.index(read[index]), 'GTAC'.index(barcode[index])
                            errors[index,read_index,barc_index] += count
                            total_count += count
                        #dists = [dist for dist in dists if dist == 1]
                        #edit_distance = sum(dist * count for dist, count in zip(dists, counts[i]))
                        #total_count += num_cycles * sum(counts[i])

                #axes[3,1].set_title('error rate {}'.format(errors / total_count * num_cycles))
                axes[3,1].plot(errors.sum(axis=(1,2)) / total_count)
                axes[3,1].set_xlabel('Cycle')
                axes[3,1].set_ylabel('Error rate')


                labels, values = np.unique(read_table.loc[~pandas.isna(read_table['aaChanges']),'edit_distance'], return_counts=True)
                assert -1 not in labels

                total_not_included = total_count - values.sum()
                labels = list(map(str, labels))
                values = list(values)
                debug ('labels = ', labels)
                debug ('values = ', values)
                #labels = ['None'] + list(map(str, labels))
                #values = [total_not_included] + list(values)
                axes[3,2].bar(labels, values, width=1, label='sing.')

                #plotting all reads that were removed in the none bar
                #axes[0,3].bar(['None'], [total_not_included], width=1)

                #axes[0,3].set_title('Edit distance, single v multiple')
                axes[3,2].set_xlabel('Combined edit distance' if double_barcode else 'Edit distance')
                axes[3,2].set_ylabel('Count')
                axes[3,2].ticklabel_format(axis='y', style='sci', scilimits=(0,0))




        nrows, ncols = 4, 4

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols, 4*nrows))

        make_plots(fig, axes)

        fig.tight_layout()
        fig.savefig(output.plot, dpi=300)

        figs = []
        axes = []
        for i in range(nrows):
            for j in range(ncols):
                subfig, subaxis = plt.subplots(figsize=(4,4))
                axes.append(subaxis)
                figs.append(subfig)

        axes = np.array(axes, dtype=object).reshape(nrows, ncols)

        make_plots(fig, axes)

        for path, fig in zip(output.plots, figs):
            fig.tight_layout()
            fig.savefig(path, dpi=300)


