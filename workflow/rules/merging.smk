import os
import glob


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


#ruleorder: merge_tables_phenotype > merge_grid_pheno_tables

def find_phenotype_table(wildcards):
    if wildcards.phenotype_type == '':
        return []
    return [phenotyping_dir + '{path_nogrid}/{phenotype_type}cells_phenotype.csv']

rule merge_phenotype_genotype:
    input:
        cell_table = segmentation_dir + '{path_nogrid}/cells.csv',
        reads = sequencing_dir + '{path_nogrid}/cells_reads.csv',
        phenotype_tables = find_phenotype_table,
    output:
        table = output_dir + '{path_nogrid}.{phenotype_type}cells_full.csv'
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 5 + 10000
    wildcard_constraints:
        path_nogrid = '((?!_grid)[^.])*',
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

