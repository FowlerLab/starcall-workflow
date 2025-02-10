import os
import sys
import glob
import time

##################################################
## Background calculation / correction
##################################################

""" So that both raw.tif and corrected.tif are in the same folder, to
make the wildcards for rules that read in both of them more simple
"""
rule link_input_stitching_raw:
    input:
        input_dir + 'well{well}/cycle{cycle}/raw.tif'
    output:
        stitching_dir + 'well{well}/cycle{cycle}/raw.tif'
    localrule: True
    shell:
        "cp -l {input[0]} {output[0]}"

""" Preforms background correction by estimating the 
"""
rule calc_background:
    input:
        lambda wildcards: expand(input_dir + 'well{well}/cycle{cycle}/raw.tif',
                well=wells, cycle=phenotype_cycles if wildcards.ispt == '_pt' else cycles),
        #expand(input_dir + 'well{well}/cycle{cycle}/raw.tif', cycle=cycles_pt, allow_missing=True),
        #input_dir + 'well{well}/cycle{cycle}/raw.tif',
    output:
        background = stitching_dir + 'background{ispt,_pt|}.tif',
        #background = stitching_dir + 'well{well}/background.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb / 2 + 5000
        #mem_mb = 1000000
    run:
        import numpy as np
        import tifffile
        #import fisseq.correction

        from pybasic.shading_correction import BaSiC

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

            #background[begin:end] = fisseq.correction.estimate_background(cur_batch_image)
            #background[begin:end] = cur_batch_image.mean(axis=0)
            background[0,begin:end] = cur_batch_image.mean(axis=0)
            background[1:,begin:end] = np.percentile(cur_batch_image, [1,5,50,95,99], axis=0)

        tifffile.imwrite(output.background, background.reshape(6, *image_shape))
        """

rule correct_background:
    input:
        images = stitching_input_dir + 'well{well}/cycle{cycle}/raw.tif',
        background = lambda wildcards: stitching_dir + 'background{}.tif'.format('_pt' if wildcards.cycle in phenotype_cycles else '')
    output:
        images = stitching_dir + 'well{well}/cycle{cycle}/corrected.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 1.5 + 10000
    run:
        import tifffile
        import fisseq.correction

        images = tifffile.imread(input.images)
        background = tifffile.imread(input.background)

        fisseq.correction.illumination_correction(images, background=background, out=images)

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
    input:
        images = stitching_input_dir + '{prefix}/cycle{cycle}/{corrected}.tif',
        positions = stitching_dir + '{prefix}/cycle{cycle}/positions.csv',
        composite = stitching_dir + '{prefix}/composite.bin',
    output:
        image = stitching_output_dir + '{prefix}/cycle{cycle}/{corrected,raw|corrected}.tif',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 2.4 + 10000
    run:
        import numpy as np
        import tifffile
        import fisseq.stitching
        import fisseq.correction
        import fisseq.utils
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
        #background = fisseq.correction.estimate_background(good_images)
        #del good_images
        #fisseq.correction.illumination_correction(images, out=images, background=background)

        full_composite = fisseq.stitching.CompositeImage.load(input.composite)
        mins = full_composite.boxes.pos1[:,:2].min(axis=0) * 2
        maxes = full_composite.boxes.pos2[:,:2].max(axis=0) * 2
        debug (mins, maxes)

        try:
            final_poses = final_poses[:,2:]

            composite = fisseq.stitching.CompositeImage()
            composite.add_images([image[:,:,0] for image in images], positions=final_poses)
            full_image = composite.stitch_images(
                real_images=images,
                mins=mins,
                maxes=maxes,
                merger=fisseq.stitching.EfficientNearestMerger(),
            )

            full_image = full_image.transpose([2,0,1])
            debug (full_image.min(), full_image.max(), full_image.dtype, full_image.shape, full_image.nbytes)
            debug (full_image.min(), full_image.max(), full_image.dtype, full_image.shape, full_image.nbytes)
            del composite
            del images
            debug (full_image.min(), full_image.max(), full_image.dtype, full_image.shape, full_image.nbytes)
            #tifffile.imwrite(output.image, full_image)
            print_mem('full_image', resources.mem_mb)
            debug (fisseq.utils.human_readable(full_image.nbytes))
            tifffile.imwrite(output.image, full_image)
        except:
            import traceback
            traceback.print_exc()
            raise

rule stitch_well:
    input:
        images = expand(stitching_input_dir + 'well{well}/cycle{cycle}/{corrected}.tif', cycle=cycles, allow_missing=True),
        positions = expand(stitching_dir + 'well{well}/cycle{cycle}/positions.csv', cycle=cycles, allow_missing=True),
        composite = stitching_dir + 'well{well}/composite.bin',
    output:
        stitching_output_dir + 'well{well}/{corrected,raw|corrected}.tif'
    resources:
        mem_mb = lambda wildcards, input: input.size_mb / 4 + 10000
    run:
        import numpy as np
        import tifffile
        import fisseq.stitching
        import fisseq.correction
        import tifffile
        import time
        import shutil

        composite = fisseq.stitching.CompositeImage.load(input.composite)
        debug(composite.boxes.pos1[:,:2].min(axis=0), composite.boxes.pos2[:,:2].max(axis=0))
        mins = composite.boxes.pos1[:,:2].min(axis=0)
        maxes = composite.boxes.pos2[:,:2].max(axis=0)
        dims = maxes - mins
        debug (mins, maxes)

        tmp_output = resources.tmpdir + '/well.tif'

        for i in fisseq.utils.simple_progress(range(len(input.images))):
            debug("stitching cycle", i, time.asctime())
            final_poses = np.loadtxt(input.positions[i], delimiter=',', dtype=int)
            debug(final_poses)

            images = tifffile.imread(input.images[i])
            images = images.transpose([0,2,3,1])

            final_poses = final_poses[:,2:]
            debug (images.shape, final_poses.shape)

            if i == 0:
                well_image = tifffile.memmap(tmp_output, shape=(len(input.images), images.shape[3], int(dims[0]), int(dims[1])), dtype=images.dtype)
            else:
                well_image = tifffile.memmap(tmp_output)
            out_full_image = well_image[i].transpose(1,2,0)

            composite = fisseq.stitching.CompositeImage()
            composite.add_images([image[:,:,0] for image in images], positions=final_poses)
            full_image = composite.stitch_images(
                real_images=images,
                mins=mins,
                maxes=maxes,
                #out=out_full_image,
                merger=fisseq.stitching.EfficientNearestMerger(),
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
    import fisseq.stitching
    import fisseq.correction
    import numpy as np

    full_image = None

    for i,(path,composite_path) in enumerate(zip(image_paths, composite_paths)):
        images = tifffile.memmap(path, mode='r')
        images = images.transpose([0,2,3,1])

        composite = fisseq.stitching.CompositeImage.load(composite_path)
        composite.images = images

        if full_image is None:
            full_image = np.zeros((len(image_paths), images.shape[3], maxes[0] - mins[0], maxes[1] - mins[1]), images.dtype)
            debug (full_image.shape)

        final_image = composite.stitch_images(
            mins = mins,
            maxes = maxes,
            merger = fisseq.stitching.EfficientNearestMerger(),
            out = full_image[i].transpose(1,2,0)
        )

    return full_image


rule stitch_well_pt:
    input:
        images_pt = expand(stitching_dir + 'well{well}/cycle{cycle}/{corrected}.tif', cycle=phenotype_cycles, allow_missing=True),
        composites_pt = expand(stitching_dir + 'well{well}/cycle{cycle}/composite.bin', cycle=phenotype_cycles, allow_missing=True),
        full_composite = stitching_dir + 'well{well}/composite.bin',
    output:
        image = stitching_output_dir + 'well{well}/{corrected,raw|corrected}_pt.tif',
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 1.5
    run:
        import numpy as np
        import tifffile
        import fisseq.stitching

        full_composite = fisseq.stitching.CompositeImage.load(input.full_composite)
        mins, maxes = full_composite.boxes.pos1.min(axis=0)[:2], full_composite.boxes.pos2.max(axis=0)[:2]
        mins *= phenotype_scale
        maxes *= phenotype_scale

        tifffile.imwrite(output.image, stitch_well_section(input.images_pt, input.composites_pt, mins, maxes))


rule stitch_well_subset:
    input:
        images = expand(stitching_dir + 'well{well}_subset{size}/cycle{cycle}/{corrected}.tif', cycle=cycles, allow_missing=True),
        composites = expand(stitching_dir + 'well{well}_subset{size}/cycle{cycle}/composite.bin', cycle=cycles, allow_missing=True),
        full_composite = stitching_dir + 'well{well}_subset{size}/composite.bin',
    output:
        image = stitching_output_dir + 'well{well}_subset{size,\d+}/{corrected,raw|corrected}.tif'
    run:
        import fisseq.stitching
        import numpy as np
        import tifffile

        full_composite = fisseq.stitching.CompositeImage.load(input.full_composite)
        mins, maxes = full_composite.boxes.pos1.min(axis=0)[:2], full_composite.boxes.pos2.max(axis=0)[:2]
        center = mins + ((maxes - mins) // 2)
        radius = int(wildcards.size) // 2
        mins, maxes = center - radius, center + radius

        tifffile.imwrite(output.image, stitch_well_section(input.images, input.composites, mins, maxes))

rule stitch_well_subset_pt:
    input:
        images = expand(stitching_dir + 'well{well}_subset{size}/cycle{cycle}/{corrected}.tif', cycle=phenotype_cycles, allow_missing=True),
        composites = expand(stitching_dir + 'well{well}_subset{size}/cycle{cycle}/composite.bin', cycle=phenotype_cycles, allow_missing=True),
        full_composite = stitching_dir + 'well{well}_subset{size}/composite.bin',
    output:
        image = stitching_output_dir + 'well{well}_subset{size,\d+}/{corrected,raw|corrected}_pt.tif'
    run:
        import fisseq.stitching
        import numpy as np
        import tifffile

        full_composite = fisseq.stitching.CompositeImage.load(input.full_composite)
        mins, maxes = full_composite.boxes.pos1.min(axis=0)[:2], full_composite.boxes.pos2.max(axis=0)[:2]
        center = mins + ((maxes - mins) // 2)
        radius = int(wildcards.size) // 2
        mins, maxes = center - radius, center + radius
        mins *= phenotype_scale
        maxes *= phenotype_scale

        tifffile.imwrite(output.image, stitch_well_section(input.images, input.composites, mins, maxes))


##################################################
## creating and stitching tiles
##################################################

rule split_grid_composite:
    input:
        composite = stitching_dir + '{prefix}/composite.bin'
    output:
        composite = stitching_output_dir + '{prefix}_seqgrid{grid_size,\d+}/grid_composite.bin',
        table = stitching_output_dir + '{prefix}_seqgrid{grid_size,\d+}/grid_positions.csv',
    resources:
        mem_mb = 5000
    run:
        import numpy as np
        import tifffile
        import fisseq.stitching

        inpcomposite = fisseq.stitching.CompositeImage.load(input.composite)
        image = np.empty((inpcomposite.boxes.pos2.max(axis=0) - inpcomposite.boxes.pos1.min(axis=0)))
        debug(image.shape)
        grid_size = int(wildcards.grid_size)

        composite = fisseq.stitching.CompositeImage()
        composite.add_split_image(image, grid_size, channel_axis=-1, overlap=cellpose_diameter*2)
        composite.save(output.composite, save_images=False)
        tile_numbers = np.arange(len(composite.boxes)).reshape(-1,1)
        np.savetxt(output.table, np.concatenate([tile_numbers // grid_size, tile_numbers % grid_size, composite.boxes.pos1], axis=1),
                delimiter=',', fmt='%d', header="tile_x,tile_y,pixel_x,pixel_y", comments='')


rule stitch_tile_well:
    input:
        images = expand(stitching_dir + '{prefix}/cycle{cycle}/{corrected}.tif', cycle=cycles, allow_missing=True),
        composites = expand(stitching_dir + '{prefix}/cycle{cycle}/composite.bin', cycle=cycles, allow_missing=True),
        grid_composite = stitching_output_dir + '{prefix}_seqgrid{grid_size}/grid_composite.bin',
    output:
        image = stitching_output_dir + '{prefix}_seqgrid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/{corrected,raw|corrected}.tif',
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 2.2 / (int(wildcards.grid_size)**2)
    run:
        import numpy as np
        import tifffile
        import fisseq.stitching

        grid_composite = fisseq.stitching.CompositeImage.load(input.grid_composite)
        grid_size, x, y = int(wildcards.grid_size), int(wildcards.x), int(wildcards.y)
        box = grid_composite.boxes[x*grid_size+y]
        debug(box)

        tifffile.imwrite(output.image, stitch_well_section(input.images, input.composites, box.pos1, box.pos2))


rule stitch_tile_well_pt:
    input:
        images_pt = expand(stitching_dir + '{prefix}/cycle{cycle}/{corrected}.tif', cycle=phenotype_cycles, allow_missing=True),
        composites_pt = expand(stitching_dir + '{prefix}/cycle{cycle}/composite.bin', cycle=phenotype_cycles, allow_missing=True),
        grid_composite = stitching_output_dir + '{prefix}_seqgrid{grid_size}/grid_composite.bin',
    output:
        image = stitching_output_dir + '{prefix}_seqgrid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/{corrected,raw|corrected}_pt.tif',
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 2.2 / (int(wildcards.grid_size)**2)
    run:
        import numpy as np
        import tifffile
        import fisseq.stitching

        grid_composite = fisseq.stitching.CompositeImage.load(input.grid_composite)
        grid_size, x, y = int(wildcards.grid_size), int(wildcards.x), int(wildcards.y)
        box = grid_composite.boxes[x*grid_size+y]
        debug(box)
        box.pos1 *= phenotype_scale
        box.pos2 *= phenotype_scale
        debug(box)

        tifffile.imwrite(output.image, stitch_well_section(input.images_pt, input.composites_pt, box.pos1, box.pos2))


