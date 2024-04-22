# FISSEQ data pipeline

This repository is a full data pipeline for the analysis of FISSEQ (Flourescent in-situ sequencing) data.
This readme will give an overview of the pipeline and how to run it, for more information there is documentation
for the fisseq python package

## Structure

The pipeline is split into two main parts, one being a python package that encapsulates the major steps
involved in FISSEQ data analysis. This is meant to be applicable to any data source, microscope, or
experiment, and tries to provide a general solution for the analysis.

The other part is a Snakemake pipeline that brings together all the functions of the python package
into a concrete data pipeline. This is specific for the analysis and experiments being done in the
Fowler lab, specifically it assumes a Nikon microscope, 8-12 sequencing cycles and 1 phenotype cycle,
4 sequencing channels, etc. If you are adapting this code for a different setup, you can modify this
pipeline or write your own using the python library.

This readme will describe running the Snakemake pipeline, for an overview of the python library see the
docs.

### General file structure:

The pipeline uses three main directories to hold the data, these being rawinput/, input/, and process/.
The specific names of these can be modified in the Snakefile if necessary.

rawinput/ is used to hold the input from the microscope, without any processing done yet.

input/ holds the image files after preprocessing and stitching, but before any FISSEQ specific
analysis has happened.

Both the input and process directories follow a similar structure, being well??.tif, well??/cycle??.tif, or well??/tile??/cycle??.tif,
depending on the step of the pipeline. Wells are the top level structure, then tiles, then cycles.

### Inputs

The input of the pipeline are the raw .nd2 files from the microscope. There should be one .nd2 file for
each well and each cycle. These files are placed in the rawinput/ directory, and the structure should
look like this:

```
rawinput/
├── 20240107_153510_433/
│   ├── Well1_ChannelDAPI,GFP,G,T,A,C_Seq0000.nd2
│   ├── Well2_ChannelDAPI,GFP,G,T,A,C_Seq0001.nd2
│   ├── Well3_ChannelDAPI,GFP,G,T,A,C_Seq0002.nd2
│   ├── Well4_ChannelDAPI,GFP,G,T,A,C_Seq0005.nd2
│   ├── Well5_ChannelDAPI,GFP,G,T,A,C_Seq0004.nd2
│   └── Well6_ChannelDAPI,GFP,G,T,A,C_Seq0003.nd2
├── 20240107_175525_477/
│   ├── Well1_ChannelDAPI,GFP,G,T,A,C_Seq0000.nd2
│   ├── Well2_ChannelDAPI,GFP,G,T,A,C_Seq0001.nd2
│   ├── Well3_ChannelDAPI,GFP,G,T,A,C_Seq0002.nd2
│   ├── Well4_ChannelDAPI,GFP,G,T,A,C_Seq0005.nd2
│   ├── Well5_ChannelDAPI,GFP,G,T,A,C_Seq0004.nd2
│   └── Well6_ChannelDAPI,GFP,G,T,A,C_Seq0003.nd2
│   ...
└── phenotype/
    ├── Well1_Channel408 nm,473 nm,545 nm,635 nm_Seq0000.nd2
    ├── Well2_Channel408 nm,473 nm,545 nm,635 nm_Seq0001.nd2
    ├── Well3_Channel408 nm,473 nm,545 nm,635 nm_Seq0002.nd2
    ├── Well4_Channel408 nm,473 nm,545 nm,635 nm_Seq0005.nd2
    ├── Well5_Channel408 nm,473 nm,545 nm,635 nm_Seq0004.nd2
    └── Well6_Channel408 nm,473 nm,545 nm,635 nm_Seq0003.nd2
```

Each folder in the rawinput directory is a cycle, the filenames don''t matter as long as they are in the
correct order alphabetically. The exception to this is the phenotype cycle, which has to be named
"phenotype". 

This is normally the same structure that the microscope saves the files as, so you can simply copy
them into the rawinput directory, and rename the phenotype cycle to "phenotype".

### Outputs

The main output of the pipeline is a csv table with sequencing information for each cell.
The structure of this file is as follows:

```
 , xpos,                ypos,               bbox_x1, bbox_y1, bbox_x2, bbox_y2, mask,
1, 4.858585858585859,   66.3989898989899,   0,       56,      12,      77,      000000000000000000001111110000000000000000000000000001111111110000...
3, 15.265079365079366,  180.40238095238095, 0,       157,     32,      207,     000000000000000000111100000000000000000000000000011111111000000000...
```

The first section contains general information about each cell, its position, bounding box, and mask.
The mask is a binary string that when reshaped into the bounding box corresponds to the cell is
a boolean mask for it.

```
count_0, barcode_0,    quality_0,    count_1, barcode_1,    quality_1,    count_2, barcode_2,    quality_2,    count_3, barcode_3,    quality_3,    count
2,       CCCACCCCCCCC, !!!r!!!!!!!!, 0,       ????????????, !!!!!!!!!!!!, 0,       ????????????, !!!!!!!!!!!!, 0,       ????????????, !!!!!!!!!!!!, 2
1,       CCCCCCACCCCC, !!!!!!n!!!!!, 1,       CCCCGCCCCGCC, !!!!l!!!!c!!, 1,       CCCCCATCCTCC, !!!w\n}!!!!!, 1,       CCCCCCCACCCA, !!!!!W!ku!!!, 4
```

The next section of the file contains the sequencing information for each cell. The top 4 reads found in each cell are shown, with 
a count for each and a quality string. The quality string is in the same format as a fastq quality string, with the actual
Phred error values estimated by how confident the base caller is in each bases correctness.

### Quickstart

To run the whole pipeline, after placing the raw .nd2 files in their locations, you can request the final files from snakemake.
If you have 3 wells, you could run the command:

```
snakemake process/well1.cells.csv process/well2.cells.csv process/well3.cells.csv
```

Most wells are to large to run at once, so you can tell snakemake to split them up by changing it to this:

```
snakemake process/well1_grid5.cells.csv process/well2_grid5.cells.csv process/well3_grid5.cells.csv
```

Where grid5 means each well will be split into a 5 by 5 grid for processing then combined back together.

If you are running the pipeline on a cluster, you can use the provided helper script run.sh, replacing snakemake with run.sh:

```
./run.sh process/well1_grid5.cells.csv process/well2_grid5.cells.csv process/well3_grid5.cells.csv
```

## Steps of the pipeline

### Extraction from the nd2 file

Inputs:

- rawinput/{date}/Well*.nd2

Outputs:

- input/well{well}/cycle{cycle}.raw.nd2
- input/well{well}/cycle{cycle}.rawpositions.csv

Extracts the images and positional information from the microscope.

### Background correction

Inputs:

- input/well{well}/cycle{cycle}.raw.nd2

Outputs:

- input/well{well}/cycle{cycle}.corrected.tif

The first step in the pipeline is background correction, where the uneven illumination of the microscope
is corrected.

### Alignment

Inputs:

- input/well{well}/cycle{cycle}.corrected.tif
- input/well{well}/cycle{cycle}.rawpositions.csv

Outputs:

- input/well{well}.composite.bin

In this step all the tiles for the well are registered together, and the global position for each
is solved. This is normally a very intensive part of the pipeline, taking up to 6 hours.
The progress can be checked in

### Stitching (full well)

Inputs:

- input/well{well}/cycle{cycle}.corrected.tif
- input/well{well}.composite.bin

Outputs:

- input/well{well}.tif

This step combines all the tiles and their global positions from registration into a single
image. This is normally not necessary as the entire well image is too large to be processed
at once, and instead the image is split into smaller tiles.

### Stitching (into tiles)

Inputs:

- input/well{well}/cycle{cycle}.corrected.tif
- input/well{well}.composite.bin

Outputs:

- input/well{well}\_grid{grid_size}/tile{x}x{y}y.tif
- (phenotype cycle) input/well{well}\_grid{grid_size}/tile{x}x{y}y/cyclePT.tif

Instead of stitching the tiles into one full image for the well, they can be stitched
into smaller tiles that are sections of the well. This is important because the whole
well can be very large and require too much memory to process at once, so we can split
it up so each tile uses a reasonable amount of memory.

It may seem strange to stitch all the tiles together only to split them back up, but there
are a couple benefits from this:

- The tiles we create here are all perfectly aligned, the tiles from the microscope are not guaranteed
  to line up, and doing this manually would require a lot of work and precision.
- We can decide how large the tiles are to maximize the memory and cpu we have access to,
  if your machine has more memory you can increase the size of the tiles.
- The overlap between tiles is very high for the ones taken on the microscope which
  wastes compute by calculating stuff for those areas multiple times. The overlap
  for these tiles is the minimum required to combine them back together.

The parameter grid_size determines how many tiles to split the well into, for example well1_grid5 means
the there will be 25 tiles in a 5 by 5 grid.

Phenotyping images are stitched separately into another file, as they are taken at a different
scale.

### Cell segmentation

Inputs:

- input/{prefix}/cyclePT.tif

Outputs:

- process/{prefix}.cells.tif
- process/{prefix}.cells_2x.tif
- process/{prefix}.nuclei.tif
- process/{prefix}.nuclei_2x.tif

One of the important processing steps is cell segmentation, which happens on the phenotyping image.
In the input and output you can see there are no real set paths, instead cell segmentation can be
run on basically any tif file that is in the input directory.

This is also one of the more intensive processing steps of the pipeline, the memory requirements
for it can get quite large. If memory becomes an issue the size of each tile can be adjusted when
splitting up the grid.

### Finding dots

Inputs:

- input/{prefix}.tif

Outputs:

- process/{prefix}.bases.csv

This is the other main processing step in the pipeline, where the flourescent dots are detected.
The output is a csv file where the first two columns are the xy position of each dot, and the rest
of the columns are the values of the dot in each of the 4 flourescent base channels.

### Calling reads

Inputs:

- process/{prefix}.cells.tif
- process/{prefix}.bases.csv

Outputs:

- process/{prefix}.cells.csv

This is the final output of the processing section of the pipeline, where the flourescent dots
are collected in cells and the sequences are read out. The structure of this file is described in the Outputs
section in more detail, but it contains the sequencing results from all the cells in the well.

### Merging grid

Inputs:

- process/well{well}_grid{grid_size}/tile{x}x{y}y.cells.csv

Outputs:

- process/well{well}_grid{grid_size}.cells.csv

If the grid was split previously, now it is combined back together to get the final cells.csv table for the
whole well. This doesn''t modify the table much, however it does update all xy positions to be global for
the well as well as adding a tile column that records which tile the cell was originally from.

