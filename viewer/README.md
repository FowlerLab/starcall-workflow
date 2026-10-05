# Dot and cell viewer

A  web app for looking at sequencing dots and cells after the Snakemake pipeline has run. It runs on the cluster and you open it in the browser on your own computer through an SSH tunnel (if you are running the pipeline elsewhere, modify `run_viewer.sh` to fit your setup to use this viewer).

Images are cropped out of the raw `.nd2` files using the solved stitching composites. Note that these images are prior to background and z-scoring. 

## Views

Pick the **well** at the top left. The list of wells follows `workflow/rules/config.smk`: `wells` from the config if it's set (a number *n* means well1…well*n*), otherwise every well found in the `rawinput/` file names.

If files needed for a well or tile are missing (for example a well that hasn't been stitched yet) they will be listed at the top of the page. Open *All expected files* to see every files the view needs, each marked present or missing. Note that these do not include the `cells_reads.csv` and `cells.csv` tables needed to index to specific reads.

**Tiled and untiled runs**: the viewer works with wells processed in grid tiles (`sequencing/<well>_grid<N>/tile02x02y/…`, the default) and with wells processed whole (targets like `sequencing/well1/cells_reads.csv`, giving `sequencing/<well>/…` and `segmentation/<well>/…`). It picks the layout for each well from the tables it finds: tiled if they exist, otherwise untiled. `--layout tiled` or `--layout untiled` forces one, for a run that has both. An untiled well has no *Tile* menu; dot indices are those of `sequencing/<well>/cells_quality.csv`, and links look like `#dot/well1/whole/123`.

**Dot view**: pick a tile and the dot's index in that tile's `cells_quality.csv`.

- *Crop*: *fit cell* (default) shows the dot's whole cell with the dot circled, or 64 px around a dot that isn't in a cell. 32, 64, 128 or 256 shows a square of that many pixels centred on the dot.
- The sequencing images form a grid: one row per cycle, one column per `sequencing_channels` entry, in order dependng on `config.yaml`.
- The phenotype images show only the `segmentation: channels` from the config, with the dot circled and the cell outline.
- The page shows the called sequence, a chastity-per-cycle chart and a z-score-per-cycle chart. The z-scores are the `values_cycleXX_[GTAC]` columns.

**Cell view**: pick a cell ID, which is the index of `cells_reads.csv` and `cells.csv`. Note that the ID alone decides the cell; the *Tile* dropdown switches to the tile the cell is actually in. *Crop* is *fit cell* (the cell plus a margin) or a square of 32–256 px centred on the cell.

- The cell's dots are marked, and a table lists its reads and dots. Click a dot to open its dot view in the table or on the dot annotation in the expanded viewing window in the left column.
- Dots come from the `dot_indicies` column of `cells_reads.csv`. Older tables don't have that column; for those, the viewer uses the `cell` column of `cells_quality.csv` instead and says so under the table.

**Both views**

- Click any image to open it in the large viewing area. Use ← / → to step through the images.

- Every image has a colour picker: sequencing channels on their column headers, phenotype channels on their captions. Images are rendered from black to that colour. Colour choices are remembered in your browser, and *Reset colours* restores the defaults.

- *Export* saves the large view, with the annotations currently shown, as one file. *PNG* is 512, 1024 or 2048 px. *SVG* embeds the image and keeps the annotations as vector lines you can edit. The file is named after what's shown (well, tile and dot index, then the image's cycle and channel). Examples: `well1_cell363284_cycle00_A.png` for cell 363284 in well1 with cycle00_A selected, or `well1_tile02x02y_dot577113_PT_DAPI.png` for dot 577113 om well1 tile02x02y with the DAPI phenotype channel selected.

- In the viewing area, the min/max sliders change the current image's contrast. *Apply to channel* uses the same range for every image of that channel (it will update it to use that range in the composite image window, and to keep it to be returned to later).

- *Contrast mode: global* uses one display range per channel for all cycles, so cycles can be compared. Calculates this by sampling from across the well (see *sample_pixels* in imaging.py for the logic used), *per image* finds the range from each crop on its own. 

- Parts of a crop that no frame imaged in a cycle are hatched grey. This happens near the edge of the well, which shifts by a few hundred pixels between cycles. Those pixesls are kept as NaN in `raw.tif`. Fully unimaged images are labelled *not imaged in cycle XX*. The header lists unimaged/partially imaged cycles for the current crop.

- *Masks and dots* sets the overlays: *None*, *Cells*, *Cells and nuclei*, *Cells and dots*, *Cells and nuclei and dots*, and each of the dot options with the dot labels, *(ids)* for the `cells_quality.csv` index or *(sequences)* for `max_seq`. The cell outline is yellow and its nucleus orange. Nuclei come from `nuclei.csv`, which is matched row for row to `cells.csv`. Their outlines are exact while `nuclei_mask_unmatched_grid<N>.tif` exists (labelled by `orig_index`), and otherwise come from the finest downscaled mask column in the table (`mask2` if present, else `mask8`). Without `nuclei.csv` and `nuclei_mask_unmatched_grid<N>.tif` the nuclei options will draw no nuclei. Similarly, the cell outlines come from `cells_mask.tif` if it exists and `cells.csv` if it does not, and will not be drawn if neither exists.

- *All cells and dots in view* (default unchecked) draws every other cell/nuclei mask and dot inside the crop, in white. The view's own cell stays yellow and its dots keep their colours. Next to the box it shows how many were found. *NOTE: Only the tables of the tile being shown are used, so cells and dots from a neighbouring tile near a tile edge aren't included.* In the large view, click a white dot to open that dot, or click inside another cell (its outline lights up on hover) to open that cell. 

- The *Composite* box stacks several images into one. Tick images in its phenotype (cycle × channel) and sequencing (cycle × channel) tables to add them to the composite. Each one is drawn from black to its channel's colour, with that channel's contrast (*Contrast mode*, or the range set with *Apply to channel*).  Colors are added together in the composite. Sequencing images are scaled to phenotype resolution when mixed with phenotype images. By default the composite is the phenotype segmentation channels (`segmentation: channels`). *Default layers* restores this composite type. Click the composite to open it in the viewing area, with the overlays, and export it. The min/max sliders are off for the composite; set each channel's range on one of its own images using the large viewing area (must hit  *Apply to channel* to make the min/max slider adjustments apply to the composite image). The file is named after its layers on export, for example `well1_cell363284_composite_PT-DAPI+PT-Ph-WGA.png`.


## Setup (once)

```bash
conda env create -f viewer/env.yaml -p ~/.conda/envs/starcall-viewer
```

If you install it somewhere else, set `VIEWER_ENV=/path/to/env`.

## Running

This site was setup to run on a grid engine server, and expects is own compute node. It uses little memory, as it reads the nd2 files from `rawinput/` when needed. Note that `/path/to/run_dir` should contain an *already completed* starcall-workflow call (at least to the stage of outputting `cells_reads.csv`, phenotyping does not have to be complete to use this site).

```bash
qlogin -l mfree=4G
cd /path/to/starcall/repo
viewer/run_viewer.sh /path/to/run_dir --port 8000
```

Note that `run_viewer.sh` will prints the node it is on and the tunnel command to run on your personal computer to view the site (fill in the user and login node sections depending on how you accessed the cluster). For example:

```
ssh -N -L 8000:localhost:8000 -J <user>@<login node> <user>@fl003.grid.gs.washington.edu
```

Then open the link it prints, `http://localhost:8000/?token=...` (see *Access* below).

If you get "port in use", chnage the port used for `--port`.

**If you can't SSH to compute nodes**, start the server with `--host 0.0.0.0`, tunnel through the login node with `ssh -N -L 8000:<node>:8000 <user>@<login node>`, and open the printed link. In this mode the token and images cross the cluster network between the login node and the compute node unencrypted, so use the default (localhost and `ssh -J`) when you can.

## Access

For security, every request needs an access token:

- The token is made on first start and saved in `~/.config/starcall-viewer/token`, readable only by you. 
- Open the printed `?token=...` link. After a restart, `http://localhost:8000` works in that browser without the token.
- To give lab members access when you're running the site, send them the link with the token privately. They also need their own tunnel to the node.
- `--token new` makes a new token, so old links and cookies stop working. `--token VALUE` uses VALUE for this run only.
- Page and dot/cell links (`#cell/well1/419534`) don't include the token. Someone without the cookie has to open the token link first, if you're sharing links in lab.

## Files it needs

These must be kept when intermediates are cleaned up:

| file | used for |
|---|---|
| `rawinput/**.nd2` | all image data |
| `config.yaml` (+ `default-config.yaml`) | channels, cycles, scales, directories |
| `stitching/<well>/cycle*/composite.json` | position of every nd2 frame in each cycle |
| `stitching/<well>/rotation.csv` | rotation applied to the frames of each cycle (`applied_deg`). Without it, no frame is rotated. This is the case for runs from before rotation correction. |
| `stitching/<well>_grid<N>/grid_composite.json` | grid tile sections |
| `sequencing/<well>_grid<N>/tile*/cells_quality.csv` | dot positions, calls, chastity, z-scores |
| `sequencing/<well>_grid<N>/tile*/cells_reads.csv` | reads and `dot_indicies` per cell |
| `segmentation/<well>_grid<N>/tile*/cells.csv` | cell bboxes and downscaled masks (`mask2` and `mask8`; the finest present is used) |
| `segmentation/<well>_grid<N>/tile*/nuclei.csv` | *optional*: nucleus outlines (`mask2` or `mask8`) |
| `segmentation/<well>_grid<N>/tile*/nuclei_mask_unmatched_grid<N>.tif` | *optional*: exact nucleus outlines |
| `segmentation/<well>_grid<N>/tile*/cells_mask.tif` | *optional*: exact outlines. Without it, outlines come from the table masks, which are blockier (2 px steps with `mask2`, 8 px with only `mask8`). |

For an untiled well the same tables are in the well's own directories, and the well composite replaces the grid composite:

| file | used for |
|---|---|
| `stitching/<well>/composite.json` | the well's bounds, which `stitch_well` crops the whole-well `raw.tif` to |
| `sequencing/<well>/cells_quality.csv`, `cells_reads.csv` | dots and reads, as above |
| `segmentation/<well>/cells.csv`, *optional* `nuclei.csv` | cells and nuclei, as above |
| `segmentation/<well>/cells_mask.tif`, `nuclei_mask_unmatched.tif` | *optional*: exact outlines. A whole-well mask is only used if it can be memory-mapped (uncompressed); otherwise outlines come from the table masks. |

The whole-well `raw.tif` / `raw_pt.tif` aren't needed in either case.

The first time a tile is viewed, its `cells_quality.csv` and `cells_reads.csv` are converted to parquet in `<run_dir>/viewer_cache/`. That takes about 10–60 s per tile, and for an untiled well as long as all its tiles together would (use `--prebuild-cache`). After that, lookups are instant.

`--prebuild-cache` converts every tile of the starting well before the server starts. `--well` sets the well the page opens on; the default is the first well with all its files present. The cache is rebuilt automatically if a table is newer than its cached copy. Delete `viewer_cache/` whenever you like.

Other options: `--cache-dir`, and `--dots-table`, `--reads-table`, `--cells-table` to point at differently named tables (for example `cells_zscored_matched_quality.csv`).

## How crops are made

- `imaging.py` reads each cycle's composite. Box *j* is nd2 frame *j*, placed by an integer (row, col) offset.
- For a crop it opens only the overlapping frames (memory-mapped with `nd2.ND2File.read_frame`, usually 1–4 frames) and pastes them.
- The paste uses the same rules as `stitch_well_section` and constitch's `EfficientNearestMerger`:
  - Tiles are clipped to the grid tile section first (for an untiled well, to the well's bounds, as in `stitch_well`).
  - Each pixel comes from the tile whose edge is farthest away, capped at 255 px.
  - Ties go to tiles in later stitching batches.
- Cycles rotated by `detect_rotation` (`applied_deg` in `rotation.csv`) have the pasted part of each frame rotated about the frame centre. This is done the same way as `starcall.rotation.rotate_frame`, with `stitching: rotation: interpolation_order`. Frames keep their shape, so their composite positions are unchanged. This needs `scipy` in the viewer env.
- Phenotype crops use the `cyclePT` composite, at `phenotype_scale / bases_scale` times the sequencing resolution.

The resulting images should match their `raw.tif` / `raw_pt.tif` counterparts exactly. If things behave abnormally, examine those images.
