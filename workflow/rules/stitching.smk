import os
import sys
import glob
import time

wildcard_constraints:
    params_alignment = params_regex('channel', 'subpix', 'solver'),

##################################################
## Background calculation / correction
##################################################

rule link_input_stitching_raw:
    """ So that both raw.tif and corrected.tif are in the same folder, to
    make the wildcards for rules that read in both of them more simple
    """
    input:
        input_dir + '{well_stitching}/cycle{cycle}/raw.tif'
    output:
        temp(stitching_dir + '{well_stitching}/cycle{cycle}/raw_tiles.tif'),
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"

rule calc_background:
    """ Preforms background correction using BaSiC:
    https://github.com/marrlab/BaSiC
    """
    input:
        lambda wildcards: expand(input_dir + '{well_stitching}/cycle{cycle}/raw.tif',
                well_stitching=wells, cycle=phenotype_cycles if wildcards.ispt == '_pt' else cycles),
        #expand(input_dir + '{well_stitching}/cycle{cycle}/raw.tif', cycle=cycles_pt, allow_missing=True),
        #input_dir + '{well_stitching}/cycle{cycle}/raw.tif',
    output:
        background = stitching_dir + 'background{ispt,_pt|}.tif',
        #background = stitching_dir + '{well_stitching}/background.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb + 5000
        #mem_mb = 1000000
    run:
        import numpy as np
        import tifffile
        #import starcall.correction

        #from pybasic.shading_correction import BaSiC
        from basicpy import BaSiC

        images = tifffile.memmap(input[0], mode='r')
        images_shape = images.shape
        images_dtype = images.dtype
        del images

        background = np.empty(images_shape[1:], np.float32)

        for chan in range(images_shape[1]):
            chan_images = np.empty((len(input), images_shape[0], *images_shape[2:]), images_dtype)

            for i, path in enumerate(input):
                images = tifffile.memmap(path, mode='r')
                chan_images[i] = images[:,chan]
                del images

            debug (chan_images.shape)
            optimizer = BaSiC(chan_images.reshape(chan_images.shape[0] * chan_images.shape[1], *chan_images.shape[2:]), verbose=True)
            optimizer.prepare()
            optimizer.run()

            background[chan] = optimizer.flatfield_fullsize.astype(np.float32)
            del chan_images

        tifffile.imwrite(output.background, background)

        """
        all_shapes = []
        dtype = None
        for path in input:
            tmp_image = tifffile.memmap(path, mode='r')#[:10]
            all_shapes.append(tmp_image.shape)
            dtype = tmp_image.dtype
            del tmp_image

        image_shape = all_shapes[0][1:]
        image_size = np.prod(image_shape)
        num_images = sum(shape[0] for shape in all_shapes)

        batch_size = max(image_size // 32, 2048)

        batch_image = np.zeros((num_images, batch_size), dtype)
        background = np.zeros((6, image_size))

        for i in range(0, image_size, batch_size):
            debug ('running batch', i, i / image_size)

            begin, end = i, min(i + batch_size, image_size)
            cur_batch_image = batch_image[:,:end-begin]
            cur_pos = 0

            for path in input:
                tmp_image = tifffile.memmap(path, mode='r')#[:10]
                debug (tmp_image.shape)
                tmp_image = tmp_image.reshape(tmp_image.shape[0], -1)
                debug (tmp_image.shape, cur_batch_image.shape, cur_batch_image[cur_pos:cur_pos+tmp_image.shape[0]].shape, tmp_image[:,begin:end].shape)
                debug (begin, end, image_size)
                cur_batch_image[cur_pos:cur_pos+tmp_image.shape[0],:end-begin] = tmp_image[:,begin:end]
                cur_pos += tmp_image.shape[0]
                del tmp_image

            #background[begin:end] = starcall.correction.estimate_background(cur_batch_image)
            #background[begin:end] = cur_batch_image.mean(axis=0)
            background[0,begin:end] = cur_batch_image.mean(axis=0)
            background[1:,begin:end] = np.percentile(cur_batch_image, [1,5,50,95,99], axis=0)

        tifffile.imwrite(output.background, background.reshape(6, *image_shape))
        """

rule correct_background:
    """ Uses the background levels estimated by BaSiC to correct the background
    of a set of image tiles
    Only used if use_corrected is true in config.yaml
    """
    input:
        images = input_dir + '{well_stitching}/cycle{cycle}/raw.tif',
        background = lambda wildcards: stitching_dir + 'background{}.tif'.format('_pt' if wildcards.cycle in phenotype_cycles else '')
    output:
        images = temp(stitching_dir + '{well_stitching}/cycle{cycle}/corrected_tiles.tif'),
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 1.5 + 10000
    run:
        import tifffile
        import starcall.correction

        images = tifffile.imread(input.images)
        background = tifffile.imread(input.background)

        starcall.correction.illumination_correction(images, background=background, out=images)

        tifffile.imwrite(output.images, images)


##################################################
## stitching tiles into full well images
##################################################

def get_background(wildcards, is_pt=False):
    if wildcards.corrected == 'corrected':
        return stitching_dir + 'background{}.tif'.format('_pt' if is_pt else '')
    return []

def get_background_pt(wildcards):
    return get_background(wildcards, True)

rule stitch_cycle:
    """ Stitches a whole well image for a single cycle. Only really useful to inspect the
    stitching, single cycle images are not used in later steps of the pipeline
    """
    input:
        images = stitching_dir + '{well_stitching}/cycle{cycle}/{corrected}_tiles.tif',
        positions = stitching_dir + '{well_stitching}/cycle{cycle}/positions.csv',
        composite = stitching_dir + '{well_stitching}/composite.json',
    output:
        image = temp(stitching_dir + '{well_stitching}/cycle{cycle}/{corrected,raw|corrected}.tif'),
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 2.4 + 10000
    run:
        import numpy as np
        import tifffile
        import constitch
        import starcall.correction
        import starcall.utils
        import tifffile

        final_poses = np.loadtxt(input.positions, delimiter=',', dtype=int)
        #final_poses = final_poses[final_poses[:,0] == int(wildcards.cycle)]
        debug(final_poses)

        """
        paths = input.images
        image = tifffile.imread(paths[0])
        #images = np.zeros((len(paths), *image.shape), dtype=image.dtype)
        images = []
        for i,path in enumerate(paths):
            images.append(tifffile.imread(path).transpose([1,2,0]))
        """

        #images = nd2.imread(input.images).transpose([0,2,3,1])
        images = tifffile.imread(input.images)
        images = images.transpose([0,2,3,1])

        debug(images.shape)
        #good_images = images[filter_edge_tiles(final_poses)]
        #background = starcall.correction.estimate_background(good_images)
        #del good_images
        #starcall.correction.illumination_correction(images, out=images, background=background)

        full_composite = constitch.load(input.composite, constraints=False)
        mins = full_composite.boxes.points1[:,:2].min(axis=0) * 2
        maxes = full_composite.boxes.points2[:,:2].max(axis=0) * 2
        debug (mins, maxes)

        try:
            final_poses = final_poses[:,2:]

            composite = constitch.CompositeImage()
            composite.add_images([image[:,:,0] for image in images], positions=final_poses)
            full_image = composite.stitch(
                real_images=images,
                mins=mins,
                maxes=maxes,
                merger=constitch.EfficientNearestMerger(),
            )

            full_image = full_image.transpose([2,0,1])
            debug (full_image.min(), full_image.max(), full_image.dtype, full_image.shape, full_image.nbytes)
            debug (full_image.min(), full_image.max(), full_image.dtype, full_image.shape, full_image.nbytes)
            del composite
            del images
            debug (full_image.min(), full_image.max(), full_image.dtype, full_image.shape, full_image.nbytes)
            #tifffile.imwrite(output.image, full_image)
            print_mem('full_image', resources.mem_mb)
            debug (starcall.utils.human_readable(full_image.nbytes))
            tifffile.imwrite(output.image, full_image)
        except:
            import traceback
            traceback.print_exc()
            raise

print (wells)

rule stitch_well:
    """ Stitches the whole well together with all sequencing cycles.
    The resulting image will have 4 dimensions: (num_cycles, num_channels, width, height)
    This image may be unreasonably large, for big wells it is recommended to use tiles
    which are stitched directly, without having to stitch the whole well (see stitch_well_tile)
    To reduce the memory usage, cycles are stitched separately and written out one by one, so
    the whole well image does not have to be held in memory.
    """
    input:
        images = expand(stitching_dir + '{well_stitching}/cycle{cycle}/{corrected}_tiles.tif', cycle=cycles, allow_missing=True),
        composite = stitching_dir + '{well_stitching}/composite{params_alignment}.json',
    output:
        temp(stitching_dir + '{well_stitching}/{corrected,raw|corrected}{params_alignment}.tif'),
    resources:
        mem_mb = lambda wildcards, input: input.size_mb / len(cycles) * 3.5 + 10000
    run:
        import numpy as np
        import tifffile
        import constitch
        import starcall.correction
        import starcall.utils
        import tifffile
        import time
        import shutil

        composite = constitch.load(input.composite, constraints=False)
        debug(composite.boxes.points1[:,:2].min(axis=0), composite.boxes.points2[:,:2].max(axis=0))
        mins = composite.boxes.points1[:,:2].min(axis=0)
        maxes = composite.boxes.points2[:,:2].max(axis=0)
        dims = maxes - mins
        debug (mins, maxes)

        tmp_output = resources.tmpdir + '/well.tif'

        for i in starcall.utils.simple_progress(range(len(input.images))):
            debug("stitching cycle", i, time.asctime())

            images = tifffile.imread(input.images[i])
            images = images.transpose([0,2,3,1])

            if i == 0:
                well_image = tifffile.memmap(tmp_output, shape=(len(input.images), images.shape[3], int(dims[0]), int(dims[1])), dtype=images.dtype)
            else:
                well_image = tifffile.memmap(tmp_output)
            out_full_image = well_image[i].transpose(1,2,0)

            subcomposite = composite.layer(i)
            full_image = subcomposite.stitch(
                real_images=images,
                mins=mins,
                maxes=maxes,
                #out=out_full_image,
                merger=constitch.EfficientNearestMerger(),
            )

            full_image = full_image.transpose([2,0,1])

            well_image[i] = full_image
            debug("written to image")
            well_image.flush()
            debug("flushed image")
            del well_image

        debug("Moving tmp file")
        shutil.move(tmp_output, output[0])
        #"""




##################################################
## stitching smaller sections, ie subsets, tiles
##################################################

def stitch_well_section(image_paths, composite_paths, mins, maxes):
    import tifffile
    import constitch
    import starcall.correction
    import numpy as np

    full_image = None

    for i,(path,composite_path) in enumerate(zip(image_paths, composite_paths)):
        images = tifffile.memmap(path, mode='r')
        images = images.transpose([0,2,3,1])

        composite = constitch.load(composite_path, constraints=False)
        composite.images = images

        if full_image is None:
            full_image = np.zeros((len(image_paths), images.shape[3], maxes[0] - mins[0], maxes[1] - mins[1]), images.dtype)
            debug (full_image.shape)

        final_image = composite.stitch(
            mins = mins,
            maxes = maxes,
            merger = constitch.EfficientNearestMerger(),
            out = full_image[i].transpose(1,2,0)
        )

    return full_image


rule stitch_well_pt:
    """ Stitches a whole well together, with all phenotyping cycles.
    Like stitch_well, the resulting image will have 4 dimensions: (num_cycles, num_channels, width, height).
    It is common for there to be only one phenotyping cycle, in which case the first dimension
    is only size 1. This is expected for the rest of the pipeline, but can cause problems if you
    try to open it with an external pipeline. To inspect individual phenotype cycles, see
    stitch_cycle.
    """
    input:
        images_pt = expand(stitching_dir + '{well_stitching}/cycle{cycle}/{corrected}_tiles.tif', cycle=phenotype_cycles, allow_missing=True),
        composites_pt = expand(stitching_dir + '{well_stitching}/cycle{cycle}/composite.json', cycle=phenotype_cycles, allow_missing=True),
        full_composite = stitching_dir + '{well_stitching}/composite.json',
    output:
        image = temp(stitching_dir + '{well_stitching}/{corrected,raw|corrected}_pt.tif'),
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 1.5
    run:
        import numpy as np
        import tifffile
        import constitch

        full_composite = constitch.load(input.full_composite, constraints=False)
        mins, maxes = full_composite.boxes.points1.min(axis=0)[:2], full_composite.boxes.points2.max(axis=0)[:2]

        # scale from base images to phenotype
        mins *= phenotype_scale
        mins //= bases_scale
        maxes *= phenotype_scale
        maxes // bases_scale

        tifffile.imwrite(output.image, stitch_well_section(input.images_pt, input.composites_pt, mins, maxes))


rule stitch_well_section:
    """ Stitches a small section of a well, with dimensions wildcards.size pixels square.
    The region is taken from the center of the well. The output image will have the
    same shape as stitch_well: (num_cycles, num_channels, width, height)
    """
    input:
        images = expand(stitching_dir + '{well_stitching}/cycle{cycle}/{corrected}_tiles.tif', cycle=cycles, allow_missing=True),
        composites = expand(stitching_dir + '{well_stitching}/cycle{cycle}/composite{params_alignment}.json', cycle=cycles, allow_missing=True),
        full_composite = stitching_dir + '{well_stitching}/composite{params_alignment}.json',
    output:
        image = temp(stitching_dir + '{well_stitching}_section{size,\d+}/{corrected,raw|corrected}{params_alignment}.tif'),
    run:
        import constitch
        import numpy as np
        import tifffile

        full_composite = constitch.load(input.full_composite, constraints=False)
        mins, maxes = full_composite.boxes.points1.min(axis=0)[:2], full_composite.boxes.points2.max(axis=0)[:2]
        center = mins + ((maxes - mins) // 2)
        radius = int(wildcards.size) // 2
        mins, maxes = center - radius, center + radius

        tifffile.imwrite(output.image, stitch_well_section(input.images, input.composites, mins, maxes))

rule stitch_well_section_pt:
    """ Stitches all phenotyping cycles for a small section of a well, with dimensions wildcards.size pixels square.
    The region is taken from the center of the well. The output image will have the
    same shape as stitch_well_pt: (num_phenotype_cycles, num_channels, width, height)
    """
    input:
        images = expand(stitching_dir + '{well_stitching}/cycle{cycle}/{corrected}_tiles.tif', cycle=phenotype_cycles, allow_missing=True),
        composites = expand(stitching_dir + '{well_stitching}/cycle{cycle}/composite.json', cycle=phenotype_cycles, allow_missing=True),
        full_composite = stitching_dir + '{well_stitching}/composite.json',
    output:
        image = temp(stitching_dir + '{well_stitching}_section{size,\d+}/{corrected,raw|corrected}_pt.tif'),
    run:
        import constitch
        import numpy as np
        import tifffile

        full_composite = constitch.load(input.full_composite, constraints=False)
        mins, maxes = full_composite.boxes.points1.min(axis=0)[:2], full_composite.boxes.points2.max(axis=0)[:2]
        center = mins + ((maxes - mins) // 2)
        radius = int(wildcards.size) // 2
        mins, maxes = center - radius, center + radius

        # scale from base images to phenotype
        mins *= phenotype_scale
        mins //= bases_scale
        maxes *= phenotype_scale
        maxes // bases_scale

        tifffile.imwrite(output.image, stitch_well_section(input.images, input.composites, mins, maxes))


##################################################
## creating and stitching tiles
##################################################

rule split_grid_composite:
    """ Creates a grid of tiles across a full well. The grid is of size
    wildcards.grid_size, with overlap specified in config.yaml as
    stitching.overlap.
    """
    input:
        composite = stitching_dir + '{well_stitching}/composite.json'
    output:
        composite = stitching_dir + '{well_stitching}_grid{grid_size,\d+}/grid_composite.json',
        table = stitching_dir + '{well_stitching}_grid{grid_size,\d+}/grid_positions.csv',
    resources:
        mem_mb = 5000
    run:
        import numpy as np
        import tifffile
        import constitch

        inpcomposite = constitch.load(input.composite, constraints=False)
        image = np.empty((inpcomposite.boxes.points1.max(axis=0) - inpcomposite.boxes.points2.min(axis=0)))
        debug(image.shape)
        grid_size = int(wildcards.grid_size)

        composite = constitch.CompositeImage()
        composite.add_split_image(image, grid_size, channel_axis=-1, overlap=cellpose_diameter*2)
        constitch.save(output.composite, composite, save_images=False)
        tile_numbers = np.arange(len(composite.boxes)).reshape(-1,1)
        np.savetxt(output.table, np.concatenate([tile_numbers // grid_size, tile_numbers % grid_size, composite.boxes.points1], axis=1),
                delimiter=',', fmt='%d', header="tile_x,tile_y,pixel_x,pixel_y", comments='')


rule stitch_tile_well:
    """ Stitches a single tile in a grid. The output image has the same
    shape as stitch_well: (num_cycles, num_channels, width, height)
    """
    input:
        images = expand(stitching_dir + '{well_stitching}/cycle{cycle}/{corrected}_tiles.tif', cycle=cycles, allow_missing=True),
        composites = expand(stitching_dir + '{well_stitching}/cycle{cycle}/composite.json', cycle=cycles, allow_missing=True),
        grid_composite = stitching_dir + '{well_stitching}_grid{grid_size}/grid_composite.json',
    output:
        image = temp(stitching_dir + '{well_stitching}_grid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/{corrected,raw|corrected}.tif'),
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 2.2 / (int(wildcards.grid_size)**2)
    run:
        import numpy as np
        import tifffile
        import constitch

        grid_composite = constitch.load(input.grid_composite, constraints=False)
        grid_size, x, y = int(wildcards.grid_size), int(wildcards.x), int(wildcards.y)
        box = grid_composite.boxes[x*grid_size+y]
        debug(box)

        tifffile.imwrite(output.image, stitch_well_section(input.images, input.composites, box.point1, box.point2))


rule stitch_tile_well_pt:
    """ Stitches a single tile in a grid. The output image has the same
    shape as stitch_well_pt: (num_phenotyping_cycles, num_channels, width, height)
    """
    input:
        images_pt = expand(stitching_dir + '{well_stitching}/cycle{cycle}/{corrected}_tiles.tif', cycle=phenotype_cycles, allow_missing=True),
        composites_pt = expand(stitching_dir + '{well_stitching}/cycle{cycle}/composite.json', cycle=phenotype_cycles, allow_missing=True),
        grid_composite = stitching_dir + '{well_stitching}_grid{grid_size}/grid_composite.json',
    output:
        image = temp(stitching_dir + '{well_stitching}_grid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/{corrected,raw|corrected}_pt.tif'),
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 2.2 / (int(wildcards.grid_size)**2)
    run:
        import numpy as np
        import tifffile
        import constitch

        grid_composite = constitch.load(input.grid_composite, constraints=False)
        grid_size, x, y = int(wildcards.grid_size), int(wildcards.x), int(wildcards.y)
        box = grid_composite.boxes[x*grid_size+y]
        debug(box)

        # scale from bases images to phenotype
        box.position *= phenotype_scale
        box.position //= bases_scale
        box.size *= phenotype_scale
        box.size //= bases_scale
        debug(box)

        tifffile.imwrite(output.image, stitch_well_section(input.images_pt, input.composites_pt, box.point1, box.point2))



##################################################
## stitching with ASHLAR
##################################################

rule stitch_well_ashlar:
    input:
        #images = lambda wildcards: [get_nd2filename(well=wildcards.well_stitching, cycle=cycle) for cycle in cycles]
        images = expand(input_dir + '{well_stitching}/cycle{cycle}/raw.tif', cycle=cycles, allow_missing=True),
        positions = expand(input_dir + '{well_stitching}/cycle{cycle}/positions.csv', cycle=cycles, allow_missing=True),
    output:
        image = stitching_dir + '{well_stitching}/raw{sigma}_ashlar.ome.tif',
    params:
        sigma = parse_param('sigma', 0),
    wildcard_constraints:
        sigma = '|_sigma\d+(.\d+)?',
    resources:
        mem_mb = 16000
    run:
        import tifffile
        import numpy as np

        outpath = output.image + '_tmp_dirs/'
        #flags = '--flip-x --flip-y -m 500'
        flags = '-m 1000'

        if params.sigma != 0:
            flags += ' --filter-sigma {}'.format(params.sigma)

        for cycle, (image_path, poses_path) in enumerate(zip(input.images, input.positions)):
            print ('writing cycle', cycle)
            images = tifffile.imread(image_path)
            poses = np.loadtxt(poses_path, delimiter=',', dtype=int)

            os.makedirs(outpath + 'cycle{:02}/'.format(cycle), exist_ok=True)

            for i in range(images.shape[0]):
                chan = 1
                path = outpath + 'cycle{:02}/chan{}_row{:03}_col{:03}.tif'.format(cycle, chan, poses[i,0], poses[i,1])
                tifffile.imwrite(path, images[i,chan])

            del images

        filepattern = ' '.join("'filepattern|{}cycle{:02}/|pattern=chan1_row{{row:03}}_col{{col:03}}.tif|overlap=0.15'".format(
                    outpath, i) for i in range(len(input.images)))
        command = '/net/fowler/vol1/home/nbradley/miniconda3/envs/stitching/bin/ashlar {} -o {} {}'.format(
                filepattern, output.image, flags)
        debug(command)

        os.system(command)

