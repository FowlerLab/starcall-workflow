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
docs at <https://fowlerlab.github.io/fisseq/fisseq.html>

### General file structure:

The pipeline uses five main directories to hold data, being input/, stitching/, sequencing/, phenotyping/, and output/.
Additionally rawinput/ is used to hold the raw images from the microscope. The specific names of each folder can
be changed in the config file. The names of the directories generally describe the purpose of each one, input/ holds
the input files to the pipeline in a microscope independent format, stitching/ holds all files related to the
stitching and alignment of the microscope images, sequencing/ contains files used to detect and call reads and
segment cells, phenotyping/ holds files that measure the visual phenotypes of cells, and output/ holds the combined
output of all these steps. The below sections go into more detail with each section of the pipeline.

All of the processing directories (stitching/, sequencing/, phenotyping/) follow a similar structure. Each has its own
input/ and output/ folder (eg stitching/output/ contains fully stitched images which are then copied into sequencing/input/
and phenotyping/input/ for further processing). In each of these folders the files contain nested directories following
a similar pattern, including well??/, well??/tile??/, and well??/tile??/cycle??/.

The first level of directories is always
organized by well, with any well specific files residing in well??/, an example would be the raw image files for
well 1, which would be at well1/raw.tif.
Well names are found from input files, either the raw .nd2 files in rawinput/ or the .tif files in input/.

The next common level of organization are tiles, which are created by dividing up a well into an arbitrary grid of smaller
sections. This is normally necessary to reduce memory requirements and allow for better parallelization, as full well images
can reach 100k pixels square. A grid of tiles is represented in files by adding on `_grid{gridsize}` to the end of the
well name, for example `well1_grid5`. However, it is necessary to mark at what point in the pipeline the grid is being used, as
sequencing and phenotyping typically use different grid sizes depending on what analysis is needed in each step. Because of
this, either `_seqgrid` or `_phenogrid` is used. Tile filenames follow the pattern `tile??x??y/` with the x and y grid position
specified. To use the previous example, the raw image file for a single tile would be `well1_seqgrid5/tile02x02y/raw.tif`.
A more in depth description of tile splitting and merging can be found at the section Tiles below.

A less common but still possible level of organization is a section of a well, which follows the pattern `well??_subset{size}`.
This is normally used to test the pipeline on a smaller part of the well. Size specifies the side length of the section in pixels,
and the section is centered in the well. If you want to get example data before committing to running the whole well, requesting
just part of the well can be useful. As before, a section of the raw images of well 1 would be at `well1_subset1000/raw.tif`.

As mentioned above, another level of organization is possible with the different cycles of imaging being separated into
folders `cycle??/`. This is not common and only occurs in the `stitching/` folder. Once stitching and alignment has been preformed
all the cycles are aligned and are stored in one multilayer tif file.

Inside any of these directories are the actual data files, which are mostly specific to each step in the pipeline. There are a few
important files that are present in all steps, such as the raw image files: `raw.tif` and `raw_pt.tif`. raw.tif is the image data
for all sequencing cycles in a single tiff file. This is a 4 dimensional file, the first being the cycle, the next being the
different channels imaged for each cycle, and the last two being the x and y dimensions of the image. `raw_pt.tif` is similar
except it contains the phenotype images. It has the same shape, with the first dimension being the different phenotype cycles,
the next being the channels imaged, and the last two being the spacial dimensions. These two images are kept separate as phenotype
cycles are usually taken at a different magnification than sequencing cycles.

Both the input and process directories follow a similar structure, being `well??.tif`, `well??/cycle??.tif`, or `well??/tile??/cycle??.tif`,
depending on the step of the pipeline. Wells are the top level structure, then tiles, then cycles.

### Nikon Microscope Inputs

If your input is from a Nikon microscope, you can copy the files directly into the rawinput directory
as specified below. If your files are in another format, skip to the General Input section to see
the format you should get your data into.

The input of the pipeline are the raw `.nd2` files from the microscope. There should be one `.nd2` file for
each well and each cycle. These files are placed in the `rawinput/` directory, and the structure should
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
them into the rawinput directory, and rename the phenotype cycle to "phenotype". If you have multiple
phenotype cycles you can name them "phenotype1", "phenotype2", or however you want as long as each begins
with "phenotype".

### General Input

If your input is not in .nd2 images, you can transform it into this format and place it in the `input/` folder.
Each well should have a directory, inside of which are directories for each cycle. Each of these subdirectories
should have a `raw.tif` file containing the raw unstitched images, and a `positions.csv` file containing the tile positions
of each of these tiles. If your microscope outputs stitched images you can place them directly in the `stitching/output/`
folder (eg `stitching/output/well1/raw.tif`), however this will only work if they are also aligned between cycles which most
microscopes will not do. It is reccomended to provide unstitched images so the stitching algorithm can align across cycles
as well as across wells at the same time.

An example `input/` folder is shown below:

```
input
├── well1
│   ├── cycle00
│   │   ├── positions.csv
│   │   └── raw.tif
│   │ ...
│   ├── cycle11
│   │   ├── positions.csv
│   │   └── raw.tif
│   └── cyclePT
│       ├── positions.csv
│       └── raw.tif
│ ...
└── well6
    ├── cycle00
    │   ├── positions.csv
    │   └── raw.tif
    │ ...
    ├── cycle11
    │   ├── positions.csv
    │   └── raw.tif
    └── cyclePT
        ├── positions.csv
        └── raw.tif
```

The shape of the `raw.tif` file should be 4 dimensional, with the first being the tiles the microscope took,
the next being the different channels each image has, and the last two being the x and y spacial dimensions.

The `positions.csv` file has 4 columns, the first two are the grid position of each tile and the next two
are the positions of each tile in pixels.

#### Auxillary Barcode Input

When performing FISSEQ experiments, it is common to use barcodes to represent a more complex change to
the cell, in the case of VISSEQ this is a certain variant. To add this information to the output
data table, you can add a `barcodes.csv` file in the folder `input/auxdata/`. All that is required is
the first column contains the barcodes that should be used to match to cells. This table will be merged
with the output table and cells that contain a barcode will have the remaining columns added to their table.

If you would like to match multiple barcodes, you can separate them with a '-' and only cells with both will
be matched.

### Outputs

The output of the pipeline is a large table containing the genotype and phenotype of all cells in the experiment.
There are generally three sections to this table, sequencing, genotype, and phenotype data.

The first couple columns contain simple cell identification information, in the format:
```
 , xpos,                ypos,               bbox_x1, bbox_y1, bbox_x2, bbox_y2
1, 4.858585858585859,   66.3989898989899,   0,       56,      12,      77
2, 15.265079365079366,  180.40238095238095, 0,       157,     32,      207
```
The position of each cell is its centroid, and the bounding box specifies the section
of image needed to contain the cell.

The first main section is the sequencing data, and contains all the reads sequenced in the cell, in the format:
```
count_0, barcode_0,    quality_0,    count_1, barcode_1,    quality_1,    count_2, barcode_2,    quality_2,    count_3, barcode_3,    quality_3,    count
2,		 TATTAATTGTGT, i_ag]_MVbV2U, 2,       TATTAATTGTTT, Offde_l`oYNA, 2,       TATTAATTGTAT, j^ge\M_fYb,_, 1,       TATTAATTGTCT, _%\hF8jg+60Z, 5
4,		 TGCTTCACTGCT, eWjngSWHYICX, 1,       TCGTTAAATTTT, eFO!a:L3V>7Q, 1,       TCGGTTACTTTT, aUIHe$Q&VJJC, 1,       TCGGTCATTTTT, d9W9`2T*`6bb, 3
```
Each read has a count, the sequence, and a quality string. Because the table has to contain a fixed number
of columns cells with less than the max number of reads will not fill all columns and have some reads
with a count of zero. The quality string is meant to approximate the quality string provided in fastq files,
ranging from '!' meaning minimal quality to '~' meaning maximum.

The next section contains any auxillary data specified in the barcodes.csv file described in the Auxillary Barcode Input
section above. An example from a visseq experiment may look like:
```
virtualBarcode,  aaChanges,  variantType,      editDistance
TATTAATTGTTT,    T224Q,      Single Missense,  0
,			,		  ,						,  -1
```
Here the virtualBarcode, aaChanges and variantType column were specified in the barcodes.csv
auxillary file, and the first column was matched to the reads. An important column is the editDistance,
which is added when merging the barcode table. Sometimes it is necessary to correct possible errors
in the recovered sequences, so the number of base changes needed to match the barcode to the read is recorded.
This makes it simple to filter on, if you only want cells that matched perfectly only select
rows with `editDistance == 0`. You will also notice that the second cell wasn't able to be matched
to any entries in the `barcode.csv` table. This is because it had multiple possible barcodes that were the
same edit distance away, so it was not possible to tell which one should be matched. When this happens
the extra columns are left blank or as NA and editDistance is set to -1.

The final section contains any phenotyping information that was calculated. This is the most varied
section of the table as it depends heavily on the experiment you are performing. If you are expecting a
simple phenotype such as intensity or cell/nuclear shape, you can use the simple feature calculation
built into the pipeline, shown below. If you are looking for a more complex phenotype or would like
to take a more unsupervised approach to finding the phenotypes you may want to use a cellprofiler
pipeline to calculate extensive features or use a vision transformer to embed each cell image. These
different methods of phenotyping cells are described more in depth in the section below on phenotyping,
and the pipeline is meant to have this part be changed as the experiment requires.

```
axis_major_length,   axis_minor_length,  cell_ch0_min,  cell_ch0_mean,  cell_ch0_max
177.12434130888718,  68.08999567385904,  1266.0,		8928.9,		   12580.0		 ...
122.72258724354198,  81.98633015396119,  1519.0,		11306.45,      15649.0
```

All of these sections are concatenated next to each other, creating the final output table. This table
is generated for each well, eg `output/well1.features.cells_full.csv`.

### Log files

Log files are kept for all jobs that create a file, and their path is the same as the file being created with `logs/` prepended. For example,
if I am trying to create the file `input/well1.composite.bin`, the log file for that would be `logs/input/well1.composite.bin.log.out` and `logs/input/well1.composite.bin.log.err`

In the case of an error snakemake will print out a message showing which job failed, as well as the log file for it.

### Quickstart

#### Installation

To get this pipeline up and running, it should be as simple as running `pip install -r requirements.txt` after cloning
this repo. Some dependencies that can cause difficulties are tensorflow which is required for cell segmentation, and
cellprofiler, which can be used for calculating phenotypes. If they are not able to be installed refer to their respective
documentations on how to best install them on your platform.

#### Running the pipeline

The pipeline uses Snakemake for organization, and the way snakemake works is you invoke it requesting certain output files.
A simple command that should work well for a standard 6 well experiment is:
```
./run.sh output/well{1..6}_seqgrid5.features.cells_full.csv
```

If you are not on a compatible SGE cluster, don't use the provided run.sh command and instead invoke snakemake directly, replacing
'./run.sh' with 'snakemake'

In the command above, there are a couple different parameters you can tweak. The first is the grid size, when I
specified seqgrid5 I indicated that the pipeline should split up the well into a 5x5 grid. You can change this to any number,
or remove it entirely to not do any splitting. You can also add phenogrid5 to split the well during phenotyping as well. For example,
if you found that memory usage was too high during phenotyping with the previous command, you could try:
```
./run.sh output/well{1..6}_seqgrid5_phenogrid5.features.cells_full.csv
```

When splitting the well into tiles, you can request only one tile:
```
./run.sh output/well{1..6}_seqgrid5/tile02x02y.features.cells_full.csv
```

Instead of only requesting one tile, you can request a small section of the well to test out the pipeline and get some example data:
```
./run.sh output/well{1..6}_subset1000.features.cells_full.csv
```

Finally, you can also change what type of phenotyping you do. If you would like to run the cellprofiler pipeline 'pipeline.cppipe'
to generate phenotype data, you can use the command:
```
./run.sh output/well{1..6}_seqgrid5_phenogrid20.cellprofiler_pipeline.cells_full.csv
```
Notice that I increased the phenotype grid to 20x20, cellprofiler does not usually run well with large images so it is necessary
to split it up more than ususal.

### Complete options for output files

As shown above, there are a lot of options that can be controlled in the file path of the output file.
I will describe them below, but it can be cumbersome to specify them this way, but by including them
in file paths we can take advantage of snakemakes dependency resolution when changing parameters. It also removes the risk of
changing parameters without rerunning the proper steps, resulting in outdated files. The different methods of
specifying parameters are described below

```
output/well{well} [_subset{size}] [_seqgrid{grid_size}] [_phenogrid{grid_size}] [/tile{x}x{y}] [.features] [.cellprofiler_{pipeline}] [.{custom_phenotyping}] .cells.csv 
```

#### Subset

When adding `_subset{size}` to a well, this means that only a size by size pixel portion of the well will be processed,
taken from the center of the well. This is useful to do a test run of the pipeline, as it greatly reduces the computation
necessary for all steps. This is similar to requesting a single tile, such as `output/well1_seqgrid20/tile09x09y/cells.csv`,
which also greatly reduces computation, but using a subset also eliminates the need to stitch the whole well together, saving
even more time.

#### Tile grids

There are two different ways to specify a grid of tiles, depending on what part of the pipeline the grid is going to
be merged at, the end of the sequencing section or the end of the phenotyping section. These parts are kept separate because merging after the
sequencing section of the pipeline is important to consolidate the different cells found during cell segmentation, making
sure there are no duplicate cells and that all cells are labeled uniquely. Grids of tiles are split creating equal size
tiles that completely cover the original image, so the resulting merged image will be exactly the same size as the original.
The overlap between tiles is possible to change, however it defaults to double the cell size, making sure that each cell
is fully in at least one tile.

The main use of tile grids is to reduce memory usage of the processing steps that will happen on the tiles. Depending on
how you are running the pipeline and what resources you have, you should adjust the grid size to avoid any memory issues.
Another benefit of splitting up the images is it allows for easy parallelization of tasks that are not multithreaded.
The main tasks like this in the pipeline currently are cell segmentation and cellprofiler, but if you add custom
phenotyping steps this may apply to them too. With these tasks many tiles can be run in parallel.

If you do not want to worry about the individual tiles you shouldn't need to, the splitting and merging are all taken
care of by the pipeline. Simply adding `_seqgrid` or `_phenogrid` into the filename will cause the input images to
be split up and merged at the end of the pipeline. However if you do want to inspect the different tiles or request
a single tile the format is quite simple, each well directory such as `well1_seqgrid5` will contain directories
`tile00x00y` up to `tile05x05y`. Each of these directories will hold the same files that a normal well directory would,
such as `raw_pt.tif`, `cells_mask.tif`, or `features.cells.csv`.

If you request a single tile output file such as `output/well1_seqgrid5_phenogrid5/tile02x02y/features.cells_full.csv`, the pipeline will
only generate the output files for that specific tile. However it will still run all tiles through the sequencing section of
the pipeline, because as explained above the sequencing grid needs to be merged back together before a phenotyping grid
can be split. If you only want to run one tile through both processing steps, you would request `output/well1_seqgrid5/tile02x02y/features.cells_full.csv`.
By not specifying a phenotype grid, the sequencing grid is never merged together and only the tile you requested is processed.
The result of this is that this resulting file can never be merged back together into a final file for the whole well because
a phenotyping grid was never created for it.

#### Phenotyping options

At the end of the filename different phenotyping methods can be selected by including certain names. These can be combined
in any way, in which case the output from each one will be concatenated together in the final table.

Including `.features` will run a simple feature calculation algorithm that returns cell shape, size area, as well as min, mean,
max, and percentile values inside the cell, nucleus, and cytoplasm for each phenotyping channel. These are enough for simple phenotype
analysis, but if you are looking for more complex phenotypes then another algorithm is probably better.

Including `.ce`



## Steps of the pipeline

### Extraction from the nd2 file

Inputs:

- `rawinput/{date}/Well*.nd2`

Outputs:

- `input/well{well}/cycle{cycle}/raw.tif`
- `input/well{well}/cycle{cycle}/positions.csv`

Extracts the images and positional information from the microscope.

### Alignment

Inputs:

- `stitching/input/well{well}/cycle{cycle}/raw.tif`
- `stitching/input/well{well}/cycle{cycle}/positions.csv`

Outputs:

- `stitching/well{well}/composite.bin`

In this step all the tiles for the well are registered together, and the global position for each
is solved. This is normally a very intensive part of the pipeline, taking up to 6 hours.
The progress can be checked in the log files for the job, at `logs/input/well{well}/composite.bin.err`

### Stitching (full well)

Inputs:

- `stitching/input/well{well}/cycle{cycle}/raw.tif`
- `stitching/well{well}/composite.bin`

Outputs:

- `stitching/output/well{well}/raw.tif`
- `stitching/output/well{well}/raw_pt.tif`

This step combines all the tiles and their global positions from registration into a single
image. This is normally not necessary as the entire well image is too large to be processed
at once, and instead the image is split into smaller tiles.

### Stitching (into tiles)

Inputs:

- `stitching/well{well}/cycle{cycle}/raw.tif`
- `stitching/well{well}/composite.bin`

Outputs:

- `stitching/output/well{well}_seqgrid{grid_size}/tile{x}x{y}y/raw.tif`
- `stitching/output/well{well}_seqgrid{grid_size}/tile{x}x{y}y/raw_pt.tif`

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

The parameter `grid_size` determines how many tiles to split the well into, for example `well1_seqgrid5` means
the there will be 25 tiles in a 5 by 5 grid.

Phenotyping images are stitched separately into another file, as they are taken at a different
scale.

### Cell segmentation

Inputs:

- `sequencing/input/{prefix}/raw_pt.tif`

Outputs:

- `sequencing/output/{prefix}/cells_mask.tif`
- `sequencing/output/{prefix}/nuclei_mask.tif`

One of the important processing steps is cell segmentation, which happens on the phenotyping image.
In the input and output you can see there are no real set paths, instead cell segmentation can be
run on basically any tif file that is in the input directory.

This is also one of the more intensive processing steps of the pipeline, the memory requirements
for it can get quite large. If memory becomes an issue the size of each tile can be adjusted when
splitting up the grid.

### Finding dots

Inputs:

- `sequencing/input/{prefix}/raw.tif`

Outputs:

- `sequencing/{prefix}/bases.csv`

This is the other main processing step in the pipeline, where the flourescent dots are detected.
The output is a csv file where the first two columns are the xy position of each dot, and the rest
of the columns are the values of the dot in each of the 4 flourescent base channels.

### Calling reads

Inputs:

- `sequencing/output/{prefix}/cells_mask.tif`
- `sequencing/{prefix}/bases.csv`

Outputs:

- `sequencing/output/{prefix}/cells.csv`

This is the final output of the processing section of the pipeline, where the flourescent dots
are collected in cells and the sequences are read out. The structure of this file is described in the Outputs
section in more detail, but it contains the sequencing results from all the cells in the well.

### Merging sequencing grid

Inputs:

- `sequencing/output/well{well}_seqgrid{grid_size}/tile{x}x{y}y/cells.csv`
- `sequencing/output/well{well}_seqgrid{grid_size}/tile{x}x{y}y/cells_mask.tif`

Outputs:

- `sequencing/output/well{well}_seqgrid{grid_size}/cells.csv`
- `sequencing/output/well{well}_seqgrid{grid_size}/cells_mask.tif`

If the grid was split previously, now it is combined back together to get the final cells.csv table for the
whole well. This doesn''t modify the table much, however it does update all xy positions to be global for
the well as well as adding a tile column that records which tile the cell was originally from.

### Splitting phenotyping grid

Inputs:

- `sequencing/output/{prefix}/cells.csv`
- `sequencing/output/{prefix}/cells_mask.tif`

Outputs:

- `phenotyping/input/{prefix}_phenogrid{grid_size}/tile{x}x{y}y/cells.csv`
- `phenotyping/input/{prefix}_phenogrid{grid_size}/tile{x}x{y}y/cells_mask.tif`

Like the sequencing step, the phenotyping step can be run on smaller tiles split from the whole well. To make this
work we have to split the files generated from the sequencing section of the pipeline, notably the cell table
and the cell mask image. It may seem redundant to join these together then split them back up, but other than
allowing for different grid sizes in these two sections there are many benefits from doing it this way. Because
the cell segmentation was run on multiple tiles with overlap, it can be difficult to reconcile the different segmentations
of the tiles in these overlapping regions. However once this is done, we can split the cells back up making sure that
each cell is in exactly one tile, with no overlap between phenotype tiles. This makes merging the results of
the phenotyping section trivial as we can just concatenate the results for each tile together, knowing that there are no
duplicate cells between tiles.

### Calculating simple features

Inputs:

- `phenotyping/input/{prefix}/raw_pt.tif`
- `phenotyping/input/{prefix}/cells.csv`
- `phenotyping/input/{prefix}/cells_mask.tif`

Outputs:

- `phenotyping/output/{prefix}/features.cells.csv`

This is the simple phenotyping solution contained in the pipeline, which calculates a handful of useful features
that can be used to determine simple phenotypes. These features include the shape, eccentricity, area and such of the
cell, and the min, mean, max, sum and different percentile values of the cell, the nucleus and the cytoplasm for each
phenotyping channel. If you are only looking at, for example intensity in channel 2 or ratio of channel 1 to channel 2 in the nucleus,
these features should be enough for you to generate useful phenotype data. However if you are looking for more complicated
phenotypes or would like to use a more unsupervised approach to find different phenotypes, this solution is probably not
enough. The main benefit of using this is the low computational cost, it runs fast and usually does not require the images
to be split into tiles, so you can run the pipeline without `_phenogrid??`.

### Calculating features with cellprofiler

Inputs:

- `phenotyping/input/{prefix}/raw_pt.tif`
- `phenotyping/input/{prefix}/cells.csv`
- `phenotyping/input/{prefix}/cells_mask.tif`

Outputs:

- `phenotyping/output/{prefix}/cellprofiler_{pipeline}.cells.csv`

Cellprofiler is a much more robust and proven method for generating many features of cells, to the level where much more complex
analysis is possible. If you would like to run a cellprofiler pipeline, you can provide one in a specific format and the 
pipeline will run it on the phenotyping images.

To get the image data into cellprofiler your pipeline should begin with a LoadData module, reading in the data from `files.csv`
located in the default input directory. This will load the different phenotype channels as CH0, CH1, CH2, ... and the cell
masks as Cells and Nuclei. The Cells mask should be converted into objects with a ConvertImageToObjects, after which any
analysis can be done using these cell objects. When exporting data with ExportToSpreadsheet, the outputs should go into
the default output folder, and it should output a file called `Cells.csv`.

Although a bit rigid, once these requirements are satesfied any typical cellprofiler analysis can be preformed, greatly
increasing the possible downstream data analysis.

### Calculating custom features


### Merging calculated features

Inputs:

- `phenotyping/output/{prefix}/cells.csv`
- `phenotyping/output/{prefix}/features.cells.csv`
- `phenotyping/output/{prefix}/cellprofiler_pipeline.cells.csv`
		(any other phenotyping tables that have been generated)

Output:

- `phenotyping/output/{prefix}/cellprofiler_pipeline.features.cells_full.csv`
															 (any additional table names would be included here)

Once all the desired phenotyping analysis has been run, the resulting features can be combined with the cell table,
creating the final table. The filename of this table reflects the different methods of feature generation used, so
by requesting a different filename you can alter which methods run.

### Merging phenotype tiles

Input:

- `phenotyping/output/{prefix}_phenogrid{grid_size}/tile{x}x{y}/features.cells_full.csv`

Output:

- `phenotyping/output/{prefix}_phenogrid{grid_size}/features.cells_full.csv`

The last step is to merge any tiles that were split for phenotyping. As described in the splitting section,
this is very simple and only consists of concatenating the different tiles and updating cell x and y positions.
If you do not care about the positional information, you can even do this step yourself with pandas.

### Output

Once all steps have finished, the output will be present in the `output/` folder, as for example
`output/well1_seqgrid5_phenogrid20.cellprofiler_pipeline.cells_full.csv`.

