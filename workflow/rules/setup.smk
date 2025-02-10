import os


##################################################
## setting up input and output dirs
##################################################

rule link_input_stitching:
    input:
        input_dir + '{file}'
    output:
        stitching_input_dir + '{file}'
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"

rule link_input_sequencing:
    input:
        stitching_output_dir + '{file}'
    output:
        sequencing_input_dir + '{file}'
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"

rule link_output_sequencing:
    input:
        sequencing_input_dir + '{prefix}/{corrected}_pt.tif'
    output:
        sequencing_output_dir + '{prefix}/{corrected,raw|corrected}_pt.tif'
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"

rule link_input_phenotyping:
    input:
        sequencing_output_dir + '{file}'
    output:
        phenotyping_input_dir + '{file}'
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"

rule link_output:
    input:
        phenotyping_output_dir + '{file}'
    output:
        output_dir + '{file}'
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"

rule link_output_nice:
    input:
        phenotyping_output_dir + '{prefix}/{type}'
    output:
        output_dir + '{prefix}.{type,[^/]+}'
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"






##################################################
## old files
##################################################

"""
rule link_old_stitching_cycle:
    input:
        'input_old/well{well}/cycle{cycle}.composite.bin',
    output:
        stitching_dir + 'well{well}/cycle{cycle}/partial_composite.bin',
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"

rule link_old_stitching_cycle_final:
    input:
        'input_old/well{well}/cycle{cycle}.fullcomposite.bin',
    output:
        stitching_dir + 'well{well}/cycle{cycle}/composite.bin',
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"

rule link_old_stitching_well:
    input:
        'input_old/well{well}.composite.bin',
    output:
        stitching_dir + 'well{well}/composite.bin',
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"

rule link_old_segmentation:
    input:
        'process_old/{prefix}_grid{gridsize}/tile{tilepos}.cells_2x.tif',
    output:
        sequencing_dir + '{prefix}_seqgrid{gridsize}/tile{tilepos}/cells_mask.tif',
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"

rule link_old_segmentation_nuclei:
    input:
        'process_old/{prefix}_grid{gridsize}/tile{tilepos}.nuclei_2x.tif',
    output:
        sequencing_dir + '{prefix}_seqgrid{gridsize}/tile{tilepos}/nuclei_mask.tif',
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"

rule link_old_dots:
    input:
        'process_old/{prefix}_grid{gridsize}/tile{tilepos}.bases.csv',
    output:
        sequencing_dir + '{prefix}_seqgrid{gridsize}/tile{tilepos}/bases.csv',
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"

rule link_old_grid_composite:
    input:
        'input_old/well{well}_grid{gridsize}/composite.bin'
    output:
        stitching_output_dir + 'well{well}_seqgrid{gridsize}/grid_composite.bin',
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"


ruleorder: link_old_stitching_cycle > align_cycle
ruleorder: link_old_stitching_cycle_final > align_well
ruleorder: link_old_stitching_well > align_well
ruleorder: link_old_segmentation > segment_cells
ruleorder: link_old_segmentation_nuclei > segment_nuclei
ruleorder: link_old_dots > find_dots
ruleorder: link_old_grid_composite > split_grid_composite

#"""

