import os
import glob
import math

##################################################
##  Extracting the images and metadata from the microscope format
##################################################

def read_nd2(path):
    import nd2
    return nd2.imread(path)

def filter_edge_tiles(positions):
    import numpy as np
    edge_tiles = []
    for x,y in positions:
        x_mask = positions[:,0] == x
        y_mask = positions[:,1] == y
        edge_tiles.append(
            np.all(positions[x_mask][:,1] >= y) or np.all(positions[x_mask][:,1] <= y)
            or np.all(positions[y_mask][:,0] >= x) or np.all(positions[y_mask][:,0] <= x)
        )

    return ~np.array(edge_tiles)

def get_nd2filename(wildcards=None, cycle=None, well=None):
    if cycle is None: cycle = wildcards.cycle
    if well is None: well = wildcards.well_base

    if cycle in phenotype_cycles:
        date = phenotype_dates[phenotype_cycles.index(cycle)]
    else:
        index = cycles.index(cycle)
        if index >= len(dates):
            raise FileNotFoundError('Error: not enough directories found in {}, needed {} to get cycle{}'.format(rawinput_dir, index+1, cycle))
        date = dates[index]

    if well[:4] == 'well':
        well = '[wW]ell' + well[4:]

    glob_path = rawinput_dir + date + '/{well}_*.nd2'.format(well=well)
    paths = glob.glob(glob_path)
    if len(paths) == 0:
        raise FileNotFoundError('Error: no nd2 files found for "{}"'.format(glob_path))
    return paths[0]


if os.path.exists(rawinput_dir):

    rule extract_nd2:
        input:
            get_nd2filename
        output:
            #expand(input_dir + '{well_base}/cycle{cycle}/tile{tile}.tif', tile=tiles, allow_missing=True)
            #input_dir + '{well_base}/cycle{cycle}/tile{tile}.tif'
            temp(input_dir + '{well_base}/tile{tile}/cycle{cycle}.tif'),
        resources:
            mem_mb = 10000
        run:
            import numpy as np
            import tifffile
            import nd2reader

            images = nd2reader.ND2Reader(input[0])
            num_channels = images.sizes['c']
            images.iter_axes = ['v','c']

            print_mem('extract_nd2', resources.mem_mb)
            #for tile in range(0,len(images)//num_channels):
            tile = int(wildcards.tile)
            image = np.array(images[tile*num_channels:(tile+1)*num_channels])
            tifffile.imwrite(output[0], image)

    rule extract_tile:
        input:
            input_dir + '{well}/cycle{cycle}.tif'
        output:
            temp(input_dir + '{well}/tile{tile}/cycle{cycle}.tif'),
        run:
            import tifffile

            images = tifffile.memmap(input[0], mode='r')
            tile = int(wildcards.tile)
            image = images[min(tile, len(images)-1)]
            tifffile.imwrite(output[0], image)

    rule extract_nd2_positions:
        """ This rule reads the nd2 files and gets the position for each tile.
        The output file has 4 columns, the first two are the grid position of the tile,
        and the second two are the estimated pixel position of each tile.
        """
        input:
            get_nd2filename
        output:
            input_dir + '{well_base}/cycle{cycle}/positions.csv'
        resources:
            mem_mb = 2000
        run:
            import numpy as np
            import tifffile
            import nd2

            with nd2.ND2File(input[0]) as images:
                full_meta = images.ome_metadata()

                if len(full_meta.images) == 1:
                    meta = full_meta.images[0].dict()
                    #debug(meta)

                    size_x = meta['pixels']['physical_size_x']# * meta['pixels']['size_x']
                    size_y = meta['pixels']['physical_size_y']# * meta['pixels']['size_y']
                    xposes = np.array([plane['position_x'] for plane in meta['pixels']['planes'] if plane['the_c'] == 0]) / size_x
                    yposes = np.array([plane['position_y'] for plane in meta['pixels']['planes'] if plane['the_c'] == 0]) / size_y
                    assert len(xposes) > 1

                else:
                    xposes, yposes = [], []

                    for i in range(len(full_meta.images)):
                        meta = full_meta.images[i].dict()
                        #debug(meta)

                        size_x = meta['pixels']['physical_size_x']# * meta['pixels']['size_x']
                        size_y = meta['pixels']['physical_size_y']# * meta['pixels']['size_y']
                        cur_xposes = [plane['position_x'] for plane in meta['pixels']['planes'] if plane['the_c'] == 0]
                        cur_yposes = [plane['position_y'] for plane in meta['pixels']['planes'] if plane['the_c'] == 0]
                        assert len(cur_xposes) == 1 and len(cur_yposes) == 1
                        xposes.append(cur_xposes[0] / size_x)
                        yposes.append(cur_yposes[0] / size_y)
                        #debug (i, cur_xposes, cur_yposes)

                    xposes, yposes = np.array(xposes), np.array(yposes)

                #positions = np.array([-xposes, -yposes]).T
                positions = np.array([yposes, -xposes]).T
                #debug (positions)
                #debug (positions.shape)
                #debug (np.array([plane['position_x'] for plane in meta['pixels']['planes']]) / size_x)

                shift_dist = max(abs(positions[0,0] - positions[1,0]), abs(positions[0,1] - positions[1,1]))
                shift_dist = np.median(np.linalg.norm(positions[1:] - positions[:-1], axis=1))
                grid_poses = np.round(positions / shift_dist).astype(int)

                grid_poses = grid_poses - grid_poses.min(axis=0).reshape(1,-1)

                #grid_poses = grid_poses - np.ceil(grid_poses.max(axis=0).reshape(1,-1) / 2)
                #debug (grid_poses.mean(axis=0), positions.mean(axis=0), (grid_poses.max(axis=0) - grid_poses.min(axis=0)).reshape(1,-1) / 2)
                debug (grid_poses.max(axis=0) / 2, grid_poses.mean(axis=0))

                positions = positions - positions.min(axis=0).reshape(1,-1)
                #positions = positions - np.round(positions.mean(axis=0)).reshape(1,-1)
                #grid_poses = grid_poses - np.round(grid_poses.mean(axis=0)).reshape(1,-1)

                debug (grid_poses, positions, grid_poses.min(axis=0), positions.min(axis=0))

                np.savetxt(output[0], np.concatenate((grid_poses, positions), axis=1), fmt='%d', delimiter=',')

    rule copy_nd2_image:
        """ Converts an nd2 file into a tiff. The resulting tiff should have 4 dimensions,
        (num_tiles, num_cycles, width, height)
        """
        input:
            get_nd2filename
        output:
            temp(input_dir + '{well_base}/cycle{cycle}/raw.tif'),
        resources:
            mem_mb = lambda wildcards, input: input.size_mb * 2 + 5000
        run:
            import tifffile
            import nd2

            image = nd2.imread(input[0])
            tifffile.imwrite(output[0], image)



##################################################
## making a smaller well input
##################################################

rule make_section:
    """ Splits an existing well in tif format into a smaller section, by taking a subset of
    tiles from the center of the well. The size parameter determines the number of tiles
    taken, a size by size grid is taken from the center of the well. This grid is
    made of tiles from the sequencing images, for phenotyping cycles the size of the
    grid is adjusted to fit the same area, as much as possible.
    This is useful to test out different alignment or stitching methods, as
    it reduces the processing needed to stitch the well dramatically.
    """
    input:
        images = expand(input_dir + '{well_nosubset}/cycle{cycle}/raw.tif', cycle=cycles_pt, allow_missing=True),
        positions = expand(input_dir + '{well_nosubset}/cycle{cycle}/positions.csv', cycle=cycles_pt, allow_missing=True),
    output:
        images = expand(input_dir + '{well_nosubset}_subset{size,\d+}/cycle{cycle}/raw.tif', cycle=cycles_pt, allow_missing=True),
        positions = expand(input_dir + '{well_nosubset}_subset{size,\d+}/cycle{cycle}/positions.csv', cycle=cycles_pt, allow_missing=True),
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 4 / len(cycles_pt) + 5000
    run:
        import tifffile
        import numpy as np

        all_poses = []
        mins, maxes = [], []
        for path in input.positions:
            poses = np.loadtxt(path, delimiter=',', dtype=int)
            cur_mins, cur_maxes = poses[:,:2].min(axis=0), poses[:,:2].max(axis=0)
            if any(path.count('cycle' + cycle) for cycle in phenotype_cycles):
                cur_mins = np.round(cur_mins * bases_scale / phenotype_scale)
                cur_maxes = np.round(cur_maxes * bases_scale / phenotype_scale)
            all_poses.append(poses)
            mins.append(cur_mins)
            maxes.append(cur_maxes)

        mins = np.min(mins, axis=0)
        maxes = np.max(maxes, axis=0)

        #images = tifffile.imread(input.images[0])
        #tile_size = int(wildcards.size) / (min(images.shape[2:]) * 0.75)
        #debug (tile_size, wildcards.size, min(images.shape[2:]))
        #tile_radius = math.ceil(tile_size / 2)
        #del images
        #tile_radius = int(wildcards.size)
        tile_size = int(wildcards.size)

        center = np.round((maxes - mins) / 2)
        #low_bound, high_bound = center - tile_radius, center + tile_radius
        low_bound = center - (tile_size // 2)
        high_bound = center + tile_size - (tile_size // 2)

        for i, poses, path in zip(range(len(all_poses)), all_poses, input.images):
            low, high = low_bound, high_bound
            if any(path.count('cycle' + cycle) for cycle in phenotype_cycles):
                low = np.round(low * phenotype_scale / bases_scale)
                high = np.round(high * phenotype_scale / bases_scale)

            images = tifffile.imread(path)
            mask = np.all((low <= poses[:,:2]) & (poses[:,:2] < high), axis=1)
            debug (path, low, high)
            debug (poses[mask,:2])

            np.savetxt(output.positions[i], poses[mask], delimiter=',', fmt='%d')
            tifffile.imwrite(output.images[i], images[mask])
            del images


rule make_noisy_well:
    """ Add gaussian noise onto an existing well.
    Creates a copy of another well, adding gaussian noise of size sigma
    to the pixel intensities of all images.
    Mostly used to stress test the stitching algorithm
    """
    input:
        images = expand(input_dir + '{well_nonoise}/cycle{cycle}/raw.tif', cycle=cycles_pt, allow_missing=True),
        positions = expand(input_dir + '{well_nonoise}/cycle{cycle}/positions.csv', cycle=cycles_pt, allow_missing=True),
    output:
        images = expand(input_dir + '{well_nonoise}_noise{size}/cycle{cycle}/raw.tif', cycle=cycles_pt, allow_missing=True),
        positions = expand(input_dir + '{well_nonoise}_noise{size}/cycle{cycle}/positions.csv', cycle=cycles_pt, allow_missing=True),
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 16 / len(cycles_pt) + 5000
    run:
        import tifffile
        import numpy as np
        import skimage.util

        variance = float(wildcards.size)

        for i, pos_path, image_path in zip(range(len(input.images)), input.positions, input.images):
            poses = np.loadtxt(pos_path, delimiter=',', dtype=int)
            images = tifffile.imread(image_path)
            images_dtype = images.dtype

            debug (np.mean(images), np.max(images), images.dtype)
            rng = np.random.default_rng(abs(hash(image_path)))

            noise = rng.normal(scale=variance, size=images.shape)
            images = images + noise
            np.clip(images, 0, np.iinfo(images_dtype).max, out=images)
            images = images.astype(images_dtype)

            debug (np.mean(images), np.max(images), images.dtype)

            #images = skimage.util.random_noise(images.astype(np.float32), mode='gaussian', seed=rng, var=variance)
            #images = images.astype(images_dtype)

            np.savetxt(output.positions[i], poses, delimiter=',', fmt='%d')
            tifffile.imwrite(output.images[i], images)

rule make_noisy_cycle_well:
    input:
        noisy_cycle = input_dir + '{well_nonoise}_noise{size}/cycle' + cycles_pt[0] + '/raw.tif',
        images = expand(input_dir + '{well_nonoise}/cycle{cycle}/raw.tif', cycle=cycles_pt[1:], allow_missing=True),
        noisy_cycle_poses = input_dir + '{well_nonoise}_noise{size}/cycle' + cycles_pt[0] + '/positions.csv',
        positions = expand(input_dir + '{well_nonoise}/cycle{cycle}/positions.csv', cycle=cycles_pt[1:], allow_missing=True),
    output:
        images = expand(input_dir + '{well_nonoise}_cyclenoise{size}/cycle{cycle}/raw.tif', cycle=cycles_pt, allow_missing=True),
        positions = expand(input_dir + '{well_nonoise}_cyclenoise{size}/cycle{cycle}/positions.csv', cycle=cycles_pt, allow_missing=True),
    resources:
        mem_mb = 1000
    run:
        for src, dest in zip(input, output):
            os.link(src, dest)

