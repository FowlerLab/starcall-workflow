# FISSEQ data pipeline: Starcall

This repository is a full data pipeline for the analysis of FISSEQ (Flourescent in-situ sequencing) data.
This readme will give an overview of the pipeline and how to run it, for more information there is documentation
for the starcall python package

## Quick start

To get the pipeline up and running quickly there is a small example dataset available here:
<https://visseq.gs.washington.edu/data_download/> under Example testing dataset.

### Installation

To install the workflow, simply clone this repository, making sure to clone recursively so that
we get both packages.

	git clone https://github.com/FowlerLab/starcall-workflow.git
	cd starcall-workflow

The packages needed for STARCall are listed in requirements.txt, and can be installed with pip.
It is recommended that you create a conda or virtual environment, as there are a good amount
of packages needed

	conda create -n ops
	conda activate ops
	# or
	# python3 -m venv ops
	# source ops/bin/activate

	pip3 install -r requirements.txt

### Download testing dataset

Once cloned, we can download the testing dataset. Although this is a very small subset of the image
data, it still will take up ~10GB of storage once extracted. In total, downloading, extracting, 
and running the pipeline will require ~20GB of data.

	wget https://visseq.gs.washington.edu/data_download/LMNA_T3_testing_image_set.tar.gz
	tar -xf LMNA_T3_testing_image_set.tar.gz
	mv LMNA_T3_testing_image_set/* ./

### Run pipeline

With the files in the correct place, we can run the pipeline with the command below.
Some steps will be memory/cpu intensive, its recommended to have 16GB of ram
when running this image set. The number of cores can be adjusted, if the process
is killed for using too much ram it may be necessary to reduce ir.
If you are on a cluster environment, the command can be modified
as shown below:

	snakemake --configfile default-config.yaml output/well1_subset3_grid.cellprofiler_022525.cells_full.csv --cores 4

For a cluster with qsub/qdel a shell script is provided

	./run.sh output/well1_subset3_grid.cellprofiler_022525.cells_full.csv --jobs 4

For slurm clusters snakemake has a built-in flag

	snakemake --slurm --configfile default-config.yaml output/well1_subset3_grid.cellprofiler_022525.cells_full.csv --cores 4

It may take a couple hours to run, depending on the machine you are running on.
Cell segmentation can especially take a long time if you don't have a gpu
available. While snakemake is running, it will print out what jobs
are currently running.

### Expected output

When it finishes, the output of the pipeline should be in output,
`output/well1_subset3_grid.cellprofiler_022525.cells_full.csv`.
An example of the
output is contained in the testing set that was downloaded as
`expected_output/well1_subset3_grid.cellprofiler_022525.cells_full.csv`.
Comparing the reads in the generated table to the reads in this output
is a good way to make sure the pipeline is running as expected.

In addition to the output the pipeline can generate summary plots, obtained getting
snakemake to generate the file `output/well1_subset3_grid/cells_reads.svg` the same
way as with the output table. These plots are also included in `expected_output/`,
and comparing these plots can make sure the pipeline ran properly. More information
on the specific plots can be found further down.

If everything ran well, you should be ready to run the pipeline with your own data.
The next section goes into more information on how the different steps work, and
how to get your data into the right format and run it.

## Structure

The pipeline is split into two main parts, one being a python package that encapsulates the major steps
involved in FISSEQ data analysis. This is meant to be applicable to any data source, microscope, or
experiment, and tries to provide a general solution for the analysis.

The other part is a Snakemake pipeline that brings together all the functions of the python package
into a concrete data pipeline. This is more specific to the analysis of VIS-seq experiments, but it
is still meant to be applicable to other in situ sequencing datasets.
If you are adapting this code for a different setup, you can modify this
pipeline or write your own using the python library.

This readme will describe running the Snakemake pipeline, for an overview of the python library see the
docs at <https://fowlerlab.github.io/starcall-docs/starcall.html>

### General file structure:

The pipeline uses six main directories to hold data, being
`input/`, `stitching/`, `segmentation/`, `sequencing/`, `phenotyping/`, and `output/`.
Additionally rawinput/ is used to hold the raw images from the microscope.
The names of the directories generally describe the purpose of each one, input/ holds
the input files to the pipeline in a microscope independent format, stitching/ holds all files related to the
stitching and alignment of the microscope images, sequencing/ contains files used to detect and call reads and
segment cells, phenotyping/ holds files that measure the visual phenotypes of cells, and output/ holds the combined
output of all these steps. The below sections go into more detail with each section of the pipeline.

All of the processing directories (`stitching/`, `segmentation/`, `sequencing/`, `phenotyping/`) follow a similar structure.
In each of these folders the files contain nested directories following
a similar pattern, including well??/, well??/tile??/, and well??/tile??/cycle??/.

The first level of directories is always
organized by well, with any well specific files residing in well??/, an example would be the raw image files for
well 1, which would be at `input/well1/raw.tif`.
Well names are found from input files, either the raw .nd2 files in rawinput/ or the .tif files in input/.

The next common level of organization are tiles, which are created by dividing up a well into an arbitrary grid of smaller
sections. This is normally necessary to reduce memory requirements and allow for better parallelization, as full well images
can reach 100k pixels square. A grid of tiles is represented in files by adding on `_grid{gridsize}` to the end of the
well name, for example `well1_grid5`. Tile filenames follow the pattern `tile??x??y/` with the x and y grid position
specified. To use the previous example, the raw image file for a single tile would be `stitching/well1_grid5/tile02x02y/raw.tif`.
A more in depth description of tile splitting and merging can be found at the section Tiles below.

A less common but still possible level of organization is a subset of a well, which follows the pattern `well??_subset{size}`.
This is normally used to test the pipeline on a smaller part of the well. Size specifies the number of microscope tiles to include,
and the section is centered in the well. If you want to get example data before committing to running the whole well, requesting
just part of the well can be useful. As before, a section of the raw images of well 1 would be at `well1_subset3/raw.tif`.

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
them into the rawinput directory, and rename the phenotype cycle to `phenotype`. If you have multiple
phenotype cycles you can name them `phenotype1`, `phenotype_20240107...`, or however you want as long as each begins
with `phenotype`.

### General Input

If your input is not in .nd2 images, you can transform it into this format and place it in the `input/` folder.
Each well should have a directory, inside of which are directories for each cycle. Each of these subdirectories
should have a `raw.tif` file containing the raw unstitched images, and a `positions.csv` file containing the tile positions
of each of these tiles. If your microscope outputs stitched images you can place them directly in the `stitching/`
folder (eg `stitching/well1/raw.tif`), however this will only work if they are also aligned between cycles which most
microscopes will not do. It is recomended to provide unstitched images so the stitching algorithm can align across cycles
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
the cell, in the case of VIS-seq this is a certain variant. To add this information to the output
data table, you can add a `barcodes.csv` file in the folder `input/auxdata/`. All that is required is
the first column contains the barcodes that should be used to match to cells. This table will be merged
with the output table and cells that contain a barcode will have the remaining columns added to their table.
An individual file for each well can also be specified by placing it in `input/well1/auxdata/`.

If you would like to match multiple barcodes, you can separate them with a '-' and only cells with both barcodes will
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
of image needed to contain the cell. All measurements here are in the scale of the phenotype images, not the
sequencing images, which is important to remember when phenotype images are taken at a different
scale. This means that cell masks can be obtained using the bbox values into the `cells_mask.tif` file and
cell images can be obtained from `raw_pt.tif` or `corrected_pt.tif`, but if retrieving images from `raw.tif` or
`corrected.tif` the bbox positions must be rescaled according to the scale of the phenotype images.

The first main section is the sequencing data, and contains all the reads sequenced in the cell, in the format:
```
num_reads, count_0, read_0,    quality_0,    count_1, read_1,    quality_1,    count_2, read_2,    quality_2,    count_3, read_3,    quality_3,    total_count
3		   2,		 TATTAATTGTGT, i_ag]_MVbV2U, 2,       TATTAATTGTTT, Offde_l`oYNA, 2,       TATTAATTGTAT, j^ge\M_fYb,_, 1,       TATTAATTGTCT, _%\hF8jg+60Z, 5
3		   4,		 TGCTTCACTGCT, eWjngSWHYICX, 1,       TCGTTAAATTTT, eFO!a:L3V>7Q, 1,       TCGGTTACTTTT, aUIHe$Q&VJJC, 1,       TCGGTCATTTTT, d9W9`2T*`6bb, 3
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
if you are trying to create the file `input/well1/composite.json`,
the log file for that would be `logs/input/well1/composite.json.err` and `logs/input/well1/composite.json.out`

In the case of an error snakemake will print out a message showing which job failed, as well as the log file for it.

### Config files

After cloning the repo and installing the pipeline, the next step to run the pipeline is
properly configuring it to work on your data. Snakemake uses yaml files for configuration,
and the file `default-config.yaml` contains all the options that the pipeline uses with documentation.
It is recommended to make a copy of this file named `config.yaml` and edit it with your
parameters.

Important parameters that should be set include:
- `phenotype_scale` and `bases_scale`, the objective used to image the phenotyping and sequencing images
- `sequencing_channels` and `phenotyping_channels`, the names of the channels imaged, it is important
to ensure the 'G', 'T', 'A', and 'C' channels are in the right order as those will be used to generate
read sequences. Additionally make note of a channel that can be used for alignment between cycles,
such as DAPI or GFP.
- `segmentation_grid_size`, `sequencing_grid_size`, and `phenotyping_grid_size` specify the size of the
grid used for the different steps of the pipeline. This depends on the size of your input images, for
a 6 well plate we found that 5 was a good grid size for segmentation and sequencing, while cellprofiler
needed a larger grid of 20.
- `stitching.channel`, this is the channel that is used for alignment, it should not change between cycles
and has to be imaged in both the sequencing and phenotyping images
- `segmentation.diameter` and `segmentation.channels` specify the inputs to Cellpose,
diameter is the estimated size of cells in pixels and channels is the nuclear and cytoplasm
channels to run cell segmentation on
- `segmentation.use_corrected` and `phenotyping.use_corrected` both determine whether
segmentation and phenotyping should be run on background corrected images, produced
by running BaSiC (<https://github.com/marrlab/BaSiC>). If your images need background correction
both of these should be set to True, but make sure to inspect the resulting images
and make sure they look good.

These highlighted parameters are the important ones to make sure are correct, but there are many
more in the config file that can be adjusted.

### Running the pipeline

The pipeline uses Snakemake for organization, and the way snakemake works is you invoke it requesting certain output files.
A simple command that should work well for a standard 6 well experiment is:
```
./run.sh output/well{1..6}_grid.features.cells_full.csv
```

If you are not on a compatible SGE cluster, don't use the provided run.sh command and instead invoke snakemake directly, replacing
'./run.sh' with 'snakemake'

In the command above, there are a couple different parameters you can tweak in the file path.
The first is `_grid`, which means that the images will be split into the grid size specified in the config file.

When splitting the well into tiles, you can request only one tile:
```
./run.sh output/well{1..6}_grid5/tile02x02y.features.cells_full.csv
```

Instead of only requesting one tile, you can request a small section of the well to test out the pipeline and get some example data:
```
./run.sh output/well{1..6}_subset3.features.cells_full.csv
```

Finally, you can also change what type of phenotyping you do. If you would like to run the cellprofiler pipeline 'pipeline.cppipe'
to generate phenotype data, you can use the command:
```
./run.sh output/well{1..6}_grid.cellprofiler_pipeline.cells_full.csv
```

### Complete options for output files

As shown above, there are a lot of options that can be controlled in the file path of the output file.
I will describe them below, but it can be cumbersome to specify them this way, but by including them
in file paths we can take advantage of snakemakes dependency resolution when changing parameters. It also removes the risk of
changing parameters without rerunning the proper steps, resulting in outdated files. The different methods of
specifying parameters are described below

```
output/well{well} [_subset{size}] [_section{size}] [_grid] [/tile{x}x{y}] [.features] [.cellprofiler_{pipeline}] .cells_full. (csv|parquet)
```

#### Subset

When adding `_subset{size}` to a well, this means that only a size by size tile section of the input microscope tiles are included,
taken from the center of the well. This is useful to do a test run of the pipeline, as it greatly reduces the computation
necessary for all steps.

#### Section

After stitching a well or subset, the final images can be cropped by adding `_section{size}` to only include
a size by size pixel section of the well, taken from the center. This can be useful if you have stitched a whole well but want to
test all downstream steps with a small section of the stitched images.

#### Tile grids

There are two different ways to specify a grid of tiles, depending on what part of the pipeline the grid is going to
be merged at, the end of the sequencing section or the end of the phenotyping section. These parts are kept separate because merging after the
sequencing section of the pipeline is important to consolidate the different cells found during cell segmentation, making
sure there are no duplicate cells and that all cells are labeled uniquely. Grids of tiles are split creating equal size
tiles that completely cover the original image, so the resulting merged image will be exactly the same size as the original.
The overlap between tiles is possible to change in the config file, as the parameter `stitching.overlap`.

The main use of tile grids is to reduce memory usage of the processing steps that will happen on the tiles. Depending on
how you are running the pipeline and what resources you have, you should adjust the grid size to avoid any memory issues.
Another benefit of splitting up the images is it allows for easy parallelization of tasks that are not multithreaded.
The main tasks like this in the pipeline currently are cell segmentation and cellprofiler, but if you add custom
phenotyping steps this may apply to them too. With these tasks many tiles can be run in parallel.

If you do not want to worry about the individual tiles you shouldn't need to, the splitting and merging are all taken
care of by the pipeline. Simply adding `_grid` into the filename will cause the input images to
be split up and merged at the end of the pipeline. However if you do want to inspect the different tiles or request
a single tile the format is quite simple, each well directory such as `well1_grid5` will contain directories
`tile00x00y` up to `tile05x05y`. Each of these directories will hold the same files that a normal well directory would,
such as `raw_pt.tif`, `cells_mask.tif`, or `features.cells.csv`.

If you request a single tile output file such as `output/well1_grid5/tile02x02y.features.cells_full.csv`, the pipeline will
only generate the output files for that specific tile. However it will still run all tiles through the sequencing section of
the pipeline, because as explained above the segmentation grid needs to be merged back together before a phenotyping grid
can be split.

#### Phenotyping options

At the end of the filename different phenotyping methods can be selected by including certain names. These can be combined
in any way, in which case the output from each one will be concatenated together in the final table.

Including `.features` will run a simple feature calculation algorithm that returns cell shape, size area, as well as min, mean,
max, and percentile values inside the cell, nucleus, and cytoplasm for each phenotyping channel. These are enough for simple phenotype
analysis, but if you are looking for more complex phenotypes then another algorithm is probably better.

Including `.cellprofiler_{pipeline}` will invoke cellprofiler with the pipeline file `{pipeline}.cppipe`, or additionally searching
for any files matching `*{pipeline}.cppipe`. A common pattern is to date different versions of pipelines, if the files
`cellprofiler_122424.cppipe` and `cellprofiler_022525.cppipe` are present and `.cellprofiler_022525` is requested,
the second pipeline will be ran. Cellprofiler is a proven method to extract many different features from cell images,
for visualization or analyis. More info on what a pipeline should look like is included in the Cellprofiler section below.

## Parameters

As explained before snakemake uses yaml files to store configuration about the pipeline and data.
However another feature that is included in the pipeline is being able to
change these parameters in the file path, by requesting different files from
snakemake. For example, the diameter that is given to cellpose is a common parameter
that needs to be tuned, and it is specified in the config file. When you generate the file
`segmentation/well1_section1000/cells_mask.tif`, it will use the diameter specified in this file.
However, if you generate the file `segmentation/well1_section1000/cells_mask_diameter50.tif`,
the diameter will be set to 50, ignoring the value in the config file. Thus you can generate many
different files with different diameter values, and compare them to find the best one.

This is true for many steps in the pipeline, and many also pass forward parameters from earlier steps.
For example, the dot detection step has three parameters, min, max, and num that specify the
range of gaussian sigmas to check when detecting sequencing dots. The final read calling step
has the parameter maxreads, and we can adjust them all by generating the file
`sequencing/well1_section1000/cells_reads_min1_max5_num10_maxreads8`.

One thing to be careful about is the order in which parameters are specified, if they
are not in the correct order they wont be split out properly. Below when the general steps
of the pipeline are explained, parameters are listed as well.

## Steps of the pipeline

### Extraction from the nd2 file

Inputs:

- `rawinput/{date}/well*.nd2`

Outputs:

- `input/well{well}/cycle{cycle}/raw.tif`
- `input/well{well}/cycle{cycle}/positions.csv`

Extracts the images and positional information from the microscope. This is only ran if the rawinput directory
is provided, if your input files are in `.tif` format they will work directly with the pipeline, placed
in the `input/` folder according to the format specified above in the Input section.

### Alignment

Inputs:

- `input/well{well}/cycle{cycle}/raw.tif`
- `input/well{well}/cycle{cycle}/positions.csv`

Outputs:

- `stitching/well{well}/composite.json`

Params:

- `channel`: The channel that is used to align between and across cycles. Can be an integer index
or one of the names of the channels specified in the config file
- `subpix`: The level to which images are aligned below 1 pixel, so if this is 8 images are aligned
to 1/8th of a pixel. The benefit of this is usually small over 16.

In this step all the tiles for the well are registered together, and the global position for each
is solved. This is normally a very intensive part of the pipeline, taking up to multiple hours.
The progress can be checked in the log files for the job, at `logs/input/well{well}/composite.json.err`

### Stitching (full well)

Inputs:

- `input/well{well}/cycle{cycle}/raw.tif`
- `stitching/well{well}/composite.json`

Outputs:

- `stitching/well{well}/raw.tif`
- `stitching/well{well}/raw_pt.tif`

This step combines all the tiles and their global positions from registration into a single
image. This is normally not necessary as the entire well image is too large to be processed
at once, and instead the image is split into smaller tiles.

### Stitching (into tiles)

Inputs:

- `input/well{well}/cycle{cycle}/raw.tif`
- `stitching/well{well}/composite.json`

Outputs:

- `stitching/well{well}_grid{grid_size}/tile{x}x{y}y/raw.tif`
- `stitching/well{well}_grid{grid_size}/tile{x}x{y}y/raw_pt.tif`

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

The parameter `grid_size` determines how many tiles to split the well into, for example `well1_grid5` means
the there will be 25 tiles in a 5 by 5 grid.

Phenotyping images are stitched separately into another file, as they are taken at a different
scale.

### Cell segmentation

Inputs:

- `stitching/{path}/raw_pt.tif`

Outputs:

- `segmentation/{path}/cells_mask.tif`
- `segmentation/{path}/nuclei_mask.tif`
- `segmentation/{path}/cells.csv`

One of the important processing steps is cell segmentation, which happens on the phenotyping image.
In the input and output you can see there are no real set paths, instead cell segmentation can be
run on basically any tif file that is in the input directory.

This is also one of the more intensive processing steps of the pipeline, the memory requirements
for it can get quite large. If memory becomes an issue the size of each tile can be adjusted when
splitting up the grid.

In addition to segmenting the cells, each cell is listed in the cells.csv table with its
x and y centroid position and its bounding box. These coordinates are in pixels in the phenotype
image.

### Merging the segmentation grid

Inputs:

- `segmentation/well{well}_grid{grid_size}/tile{x}x{y}y/cells_mask.tif`
- `segmentation/well{well}_grid{grid_size}/tile{x}x{y}y/nuclei_mask.tif`
- `segmentation/well{well}_grid{grid_size}/tile{x}x{y}y/cells.csv`

Outputs:

- `segmentation/well{well}_grid/cells_mask.tif`
- `segmentation/well{well}_grid/nuclei_mask.tif`
- `segmentation/well{well}_grid/cells.csv`

The segmentation is merged into a full well image, combining any cells that were in the
overlapping regions between tiles.

### Splitting the segmentation grid

Inputs:

- `segmentation/{path}_grid/cells.csv`
- `segmentation/{path}_grid/cells_mask.tif`

Outputs:

- `segmentation/{path}_grid{grid_size}/tile{x}x{y}y/cells.csv`
- `segmentation/{path}_grid{grid_size}/tile{x}x{y}y/cells_mask.tif`

Like the segmentation step, the sequencing and phenotyping steps can be run on smaller tiles split from the whole well. To make this
work we have to split the files generated from the segmentation section of the pipeline.
It may seem redundant to join these together then split them back up, but other than
allowing for different grid sizes in these two sections there are many benefits from doing it this way. Because
the cell segmentation was run on multiple tiles with overlap, it can be difficult to reconcile the different segmentations
of the tiles in these overlapping regions. However once this is done, we can split the cells back up making sure that
each cell is in exactly one tile, with no overlap between tiles. This makes merging the results of
the sequencing and phenotyping sections trivial as we can just concatenate the results for each tile together,
knowing that there are no duplicate cells between tiles.

### Calling reads

Inputs:

- `segmentation/{path}/cells_mask.tif`
- `stitching/{path}/raw.tif`

Outputs:

- `sequencing/{path}/cells_raw_reads.csv`

This is the other main processing step in the pipeline, where the flourescent dots are detected,
and their sequencing values are read out. This table contains each sequencing colony that was
detected, with its position, raw image values, sequence, and cell that it is contained in.

### Combining reads

Inputs:

- `segmentation/{path}/cells_mask.tif`
- `sequencing/{path}/cells_raw_reads.csv`

Outputs:

- `sequencing/{path}/cells_reads_partial.csv`

Typically there are many of the same reads next to each other, and to combine these reads
we cluster reads together. Right now we only cluster reads that have the same sequence
and are in the same cell, but by changing the config params in the section `read_clustering`,
reads can be combined in many different ways. Once reads with the same sequence have been
combined, we aggregate all the reads in each cell together. The resulting table has a row
for each cell, and has the top 5 reads found in the cell.

### Matching reads to barcodes

Inputs:

- `input/auxdata/barcodes.csv`
- `sequencing/{path}/cells_reads_partial.csv`

Outputs:

- `sequencing/{path}/cells_reads.csv`

Here we link reads in cells to the barcode lookup table provided. For each cell, we search
for the match between one of the reads and a barcode in the lookup table with minimal edit distance.
If there is a single match with minimum edit distance, we link that cell with the barcode.
If there are multiple matches with the same edit distance, it is ambiguous which barcode
the cell should be linked to. This can happen for multiple reasons, errors in sequencing may
have caused the barcode to change enough where it is now equal edit distance between two barcodes.
Alternatively the cell could have two reads that both map to a barcode, either because
multiple barcodes are actually present in the cell, or because the segmentation of the cell
is not perfectly accurate and a neighboring cells reads are being misassigned. Either way,
we are unsure which barcode to map to and we dont link the cell with either.

Some experiments actually expect there to be multiple barcodes in cells, and use them
to map cells more accurately with less cycles. In this case, the barcodes.csv file should
have the sets of barcodes separated by '-'. Matches between sets of reads and sets of barcodes
are searched for, with the individual edit distances added together. The same procedure
is used, where if multiple matches have the same edit distance, the cell cannot be linked.

The result of this process is that some cells are linked to a barcode and thus a row in the
barcode table. All the columns of this barcode table are added to the reads table, and the values
for linked barcodes are copied to their respective cells. 

### Merging phenotype tables

Input:

- `sequencing/{path}_grid{grid_size}/tile{x}x{y}/cells_reads.csv`

Output:

- `phenotyping/{path}_grid{grid_size}/cells_reads.csv`

The last step is to merge any tiles that were split for sequencing. As described in the splitting section,
this is very simple and only consists of concatenating the different tables together. Because
no cells are contained in multiple tiles, we know there will be no duplicates.

### Calculating simple features

Inputs:

- `phenotyping/input/{path}/raw_pt.tif`
- `phenotyping/input/{path}/cells.csv`
- `phenotyping/input/{path}/cells_mask.tif`

Outputs:

- `phenotyping/{path}/features.cells.csv`

This is the simple phenotyping solution contained in the pipeline, which calculates a handful of useful features
that can be used to determine simple phenotypes. These features include the shape, eccentricity, area and such of the
cell, and the min, mean, max, sum and different percentile values of the cell, the nucleus and the cytoplasm for each
phenotyping channel. If you are only looking at, for example intensity in channel 2 or ratio of channel 1 to channel 2 in the nucleus,
these features should be enough for you to generate useful phenotype data. However if you are looking for more complicated
phenotypes or would like to use a more unsupervised approach to find different phenotypes, this solution is probably not
enough. The main benefit of using this is the low computational cost, it runs fast and usually does not require the images
to be split into tiles, so you can set the `phenotyping_grid_size` parameter in the config file to 1

### Calculating features with cellprofiler

Inputs:

- `phenotyping/input/{path}/raw_pt.tif`
- `phenotyping/input/{path}/cells.csv`
- `phenotyping/input/{path}/cells_mask.tif`

Outputs:

- `phenotyping/{path}/cellprofiler_{pipeline}.cells.csv`

Cellprofiler is a much more robust and proven method for generating many features of cells, to the level where much more complex
analysis is possible. If you would like to run a cellprofiler pipeline, you can provide one in a specific format and the 
pipeline will run it on the phenotyping images.

To get the image data into cellprofiler your pipeline should begin with a LoadData module, reading in the data from `files.csv`
located in the default input directory. This will load the different phenotype channels as CH0, CH1, CH2, ... and the cell
masks as Cells and Nuclei. The Cells mask should be converted into objects with a ConvertImageToObjects, after which any
analysis can be done using these cell objects. When exporting data with ExportToSpreadsheet, the outputs should go into
the default output folder, and it should output a file called `Cells.csv`.

Although a bit rigid, once these requirements are satisfied any typical cellprofiler analysis can be preformed, greatly
increasing the possible downstream data analysis.

### Merging calculated features

Inputs:

- `phenotyping/{path}/cells.csv`
- `phenotyping/{path}/features.cells.csv`
- `phenotyping/{path}/cellprofiler_pipeline.cells.csv`
		(any other phenotyping tables that have been generated)

Output:

- `phenotyping/{path}/cellprofiler_pipeline.features.cells_phenotype.csv`
															 (any additional table names would be included here)

Once all the desired phenotyping analysis has been run, the resulting features can be combined with the cell table,
creating the final table. The filename of this table reflects the different methods of feature generation used, so
by requesting a different filename you can alter which methods run.

### Merging phenotype tables

Input:

- `phenotyping/{path}_grid{grid_size}/tile{x}x{y}/features.cells_phenotype.csv`

Output:

- `phenotyping/{path}_grid{grid_size}/features.cells_phenotype.csv`

The last step is to merge any tiles that were split for phenotyping. As described in the splitting section,
this is very simple and only consists of concatenating the different tables together. Because
no cells are contained in multiple tiles, we know there will be no duplicates.

### Merging phenotype and genotype

Input:

- `segmentation/{path}/cells.csv`
- `sequencing/{path}/cells_reads.csv`
- `phenotyping/{path}/features.cells_phenotype.csv`

Output:

- `output/{path}.features.cells_full.csv`

The final step is merging the phenotype and genotype of cells into the final table. This is simply a join on
the cell ids shared between the tables. Every cell in the cells.csv file will have a row in the final table,
if the cell has no reads or phenotyping those columns will be missing.

## Quality control steps

There are multiple ways to make sure that the different steps in the pipeline are performing as desired.
Many of these are additional processing steps in the pipeline that can be run by requesting
the specific file with snakemake. This combined with the ability to change parameters in file paths mean
that finding the right values for different parameters is quite simple. For example, if you want to try a couple
different possible values for the diameter given to cellpose, you could request the followin files as so:

	snakemake output/qc/well1_section1000/cells_overlay_diameter{25,50,75,100}.png

Upon inspecting these 4 files, you can select the best diameter value and set it in the config file,
making it the default. This can be done with many of the qc plots and images below.

### Stitching plots

Outputs:

- `output/qc/{path}/cycle{cycle}_cycle{cycle}_scores_calculated.png`
- `output/qc/{path}/cycle{cycle}_cycle{cycle}_scores_filtered.png`
- `output/qc/{path}/presolve.png`
- `output/qc/{path}/solve.png`
- `output/qc/{path}/solve_accuracy.png`

Params:
- `channel`: The channel that is used to align between and across cycles. Can be an integer index
or one of the names of the channels specified in the config file
- `subpix`: The level to which images are aligned below 1 pixel, so if this is 8 images are aligned
to 1/8th of a pixel. The benefit of this is usually small over 16.

These plots are generated while stitching happens, and show the pairwise alignment between
overlapping tiles. The stitching algorithm works by calculating many of these pairwise alignments
between neighboring tiles in the same cycle and overlapping tiles in different cycles.
They are plotted as an arrow placed between the two tiles, pointing in the direction of the
estimated alignment, and marked with a dot that is colored by its score.

The plots between cycles show the alignments of only those two cycles, and the presolve plot
shows all alignments before they are all solved. When they are solved, the stitching algorithm
finds global positions for each image, and these are shown in the solve.png plot. solve_accuracy.png
shows the same image alignments, but they are colored instead by the error in the estimated alignment.

When looking at these plots, it can be quite hard to tell if the stitching worked well. The best plot
to inspect is the solve_accuracy.png plot, making sure that the images are in the general shape of a well.
If there is a catastrophic failure, it will be aparent, and the images will be scattered all over or
bunched up in one location.

### Read quality plots

Inputs:

- `sequencing/{path}/cells_reads.csv`

Output:

- `output/{path}/cells_reads.svg`

Params:

- min, max, num: The range of sigmas to search for blobs, passed to `skimage.feature.blob_log`
- norm, posweight, valweight, seqweight: Parameters used to calculate distances between reads, used for clustering.
- linkage, thresh: Parameters used to combine reads using the distances calculated.
- max_reads: The number of reads to keep per cell.

This figure contains a collection of quality control plots based on the reads generated.
The first row shows general stats including nucleotide frequencies, read count per cell, and edit distances
of matches. The next row shows stats on error rates, based on the mappings made to the library.
Note that these error rates are not true error rates, as they are only the mappings that were
able to be made to the library with minimal edit distance.
The next row contains a very useful plot which shows the edit distance to the best matching
barcode and the second best matching barcode. This plot shows what percent of cells didnt
have an ambiguous match, as all cells that are not on the diagonal were able to be matched.
The final row contains a couple more miscelanious plots showing edit distance and total error
rate over cycles.

### Cell segmentation overlay

Input:
	
- `stitching/{path}/raw_pt.tif`
- `segmentation/{path}/cells_mask.tif`

Output:

- `output/qc/{path}/cells_overlay.tif`
- `output/qc/{path}/cells_overlay.png`

Params:

- diameter: approximate pixel size of the cells
- nuclearchannel: the channel to use as the nuclear stain
- cytochannel: the channel to use as the cytoplasm stain

This generates an overlay image with the phenotype image that was used for cell segmentation and
the resulting labels. This can be useful to ensure cell segmentation is working, as all steps rely
on high quality cell segmentation.

### Highlight detected sequencing colonies

Input:

- `stitching/{path}/raw.tif`
- `sequencing/{path}/bases.csv`

Output:

- `output/qc/{path}/annotated.tif`

Params:
	min, max, num: The range of sigmas to search for blobs, passed to skimage.feature.blob_log

This labels each dot detected in the sequencing images by drawing an x at its location in an extra channel.
This is very useful to make sure that dots are being identified properly.

