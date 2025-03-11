import os
import glob


def get_grid_filenames_pheno(wildcards):
    grid_size = int(wildcards.grid_size)
    numbers = ['{:02}'.format(i) for i in range(grid_size)]
    return expand(phenotyping_dir + '{possible_output_dir}{prefix}_phenogrid{grid_size}/tile{x}x{y}y/{type}.csv', x=numbers, y=numbers, allow_missing=True)

def get_composite_pheno(wildcards):
    prefix = wildcards.prefix.split('_seqgrid')[0] + '_seqgrid{grid_size}'
    return stitching_output_dir + prefix + '/grid_composite.json'

def merge_csv_files(input_paths, output_path, extra_columns=tuple(), row_func=None):
    import csv

    readers = [csv.DictReader(open(path, newline='')) for path in input_paths]
    fieldnames = []
    for reader in readers:
        for name in reader.fieldnames:
            if name not in fieldnames: fieldnames.append(name)

    for name in ['file_index', 'file_path', *extra_columns]:
        if name not in fieldnames: fieldnames.append(name)

    with open(output_path, 'w', newline='') as ofile:
        writer = csv.DictWriter(ofile, fieldnames)
        writer.writeheader()

        for i, reader, path in zip(range(len(readers)), readers, input_paths):
            for row in reader:
                row['file_index'] = i
                row['file_path'] = path
                for name in extra_columns: row[name] = i
                if row_func is not None:
                    result = row_func(row)
                    if result is not None:
                        row.update(result)
                writer.writerow(row)


def get_other_tables(wildcards):
    return [phenotyping_output_dir + '{prefix}/' + name + '.csv' for name in wildcards.phenotype_tables.split('.')]

rule merge_tables_phenotype:
    input:
        cell_table = phenotyping_input_dir + '{prefix}/cells.csv',
        other_tables = get_other_tables,
    output:
        table = phenotyping_output_dir + '{prefix}/{phenotype_tables}.cells_phenotype.csv',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 2 + 10000
    run:
        import pandas

        cell_table = pandas.read_csv(input.cell_table, index_col=0)
        index = cell_table.index
        #cell_table = cell_table.reset_index(drop=True)
        #table = pandas.concat([cell_table] + [pandas.read_csv(path).iloc[:len(cell_table.index),:] for path in input.other_tables], axis=1)
        table = pandas.concat([pandas.read_csv(path).iloc[:len(cell_table.index),:] for path in input.other_tables], axis=1)
        table = table.set_index(index)
        table.to_csv(output.table)


rule merge_grid_pheno_tables:
    input:
        tables = get_grid_filenames_pheno,
        composite = get_composite_pheno,
    output:
        table = phenotyping_dir + '{possible_output_dir}{prefix}_phenogrid{grid_size,\d+}/{type,[^/]*}.csv',
    resources:
        #mem_mb = lambda wildcards, input: input.size_mb * 50 + 5000
        mem_mb = 5000
    wildcard_constraints:
        possible_output_dir = '(' + output_dir + ')|',
    run:
        import constitch

        composite = constitch.load(input.composite)

        def row_func(row):
            i = row['file_index']
            box = composite.boxes[i]
            row['tile_index'] = i
            row['tile_x'] = i // int(wildcards.grid_size)
            row['tile_y'] = i % int(wildcards.grid_size)
            row['xpos'] = float(row['xpos']) + box.position[0]
            row['ypos'] = float(row['ypos']) + box.position[1]
            row['bbox_x1'] = int(row['bbox_x1']) + box.position[0]
            row['bbox_x2'] = int(row['bbox_x2']) + box.position[0]
            row['bbox_y1'] = int(row['bbox_y1']) + box.position[1]
            row['bbox_y2'] = int(row['bbox_y2']) + box.position[1]

        merge_csv_files(input.tables, output.table, extra_columns=['tile_x', 'tile_y', 'tile_index'], row_func=row_func)

#ruleorder: merge_tables_phenotype > merge_grid_pheno_tables

def find_phenotype_table(wildcards):
    if wildcards.phenotype_type == '':
        return []
    return [phenotyping_output_dir + '{prefix}{possible_phenogrid}/{phenotype_type}cells_phenotype.csv']

rule merge_phenotype_genotype:
    input:
        cell_table = sequencing_output_dir + '{prefix}/cells.csv',
        reads = sequencing_output_dir + '{prefix}/cells_reads.csv',
        phenotype_tables = find_phenotype_table,
    output:
        table = output_dir + '{prefix}{possible_phenogrid}.{phenotype_type}cells_full.csv'
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 2 + 10000
    wildcard_constraints:
        prefix = '((?!_phenogrid).)*',
        possible_phenogrid = '(_phenogrid\d+)|',
        phenotype_type = '[^/]*',
    run:
        import pandas

        #tables = [pandas.read_csv(path, index_col=0) for path in input]
        #table = pandas.concat(tables, axis=1)
        table = pandas.read_csv(input[0], index_col=0)
        for path in input[1:]:
            table = table.join(pandas.read_csv(path, index_col=0))
        table.to_csv(output.table)

rule merge_all_well_tables:
    input:
        #tables = expand(output_dir + 'well{well}{path}.csv', well=wells, allow_missing=True),
        tables = lambda wildcards: expand(output_dir + 'well{well}{path}.csv', well=wildcards.wells.split('-'), allow_missing=True),
    output:
        table = output_dir + 'well{wells}{path}.csv',
    resources:
        mem_mb = 5000
    wildcard_constraints:
        # sequence of well names joined by '-'
        wells = '({well_pat})(-({well_pat}))+'.format(well_pat='|'.join(wells)),
    run:
        merge_csv_files(input.tables, output.table, extra_columns=['well'],
                row_func=lambda row: {'well': wells[row['file_index']]})

rule convert_to_parquet:
    input:
        table = output_dir + '{path}.csv',
    output:
        table = output_dir + '{path}.parquet',
    resources:
        mem_mb = 50000
    run:
        import pandas

        chunksize = 100000
        with pandas.read_csv(input.table, chunksize=chunksize) as table_chunks:
            for i, chunk in enumerate(table_chunks):
                print('Chunk', i)

                if i==0:
                    chunk.to_parquet(output.table,
                                engine='fastparquet',
                                row_group_offsets=chunksize,
                                index=False)
                else:
                    chunk.to_parquet(output.table,
                                engine='fastparquet',
                                row_group_offsets=chunksize,
                                index=False, 
                                append=True)
                del chunk

