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

def get_nd2filename(wildcards):
    if wildcards.cycle in phenotype_cycles:
        date = phenotype_dates[phenotype_cycles.index(wildcards.cycle)]
    else:
        index = cycles.index(wildcards.cycle)
        if index >= len(dates):
            return ['not_found']
        date = dates[index]
    path = rawinput_dir + date + '/Well{well}_*.nd2'.format(**wildcards)
    paths = glob.glob(path)
    if len(paths) == 0:
        return ['not_found']
    return paths[0]


rule extract_nd2:
    input:
        get_nd2filename
    output:
        #expand(input_dir + 'well{well}/cycle{cycle}/tile{tile}.tif', tile=tiles, allow_missing=True)
        #input_dir + 'well{well}/cycle{cycle}/tile{tile}.tif'
        input_dir + 'well{well}/tile{tile}/cycle{cycle}.tif',
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
        input_dir + 'well{well}/cycle{cycle}.tif'
    output:
        input_dir + 'well{well}/tile{tile}/cycle{cycle}.tif'
    run:
        import tifffile

        images = tifffile.memmap(input[0], mode='r')
        tile = int(wildcards.tile)
        image = images[min(tile, len(images)-1)]
        tifffile.imwrite(output[0], image)

""" This rule reads the nd2 files and gets the position for each tile.
The output file has 4 columns, the first two are the grid position of the tile,
and the second two are the estimated pixel position of each tile.
"""
rule extract_nd2_positions:
    input:
        get_nd2filename
    output:
        input_dir + 'well{well}/cycle{cycle}/positions.csv'
    resources:
        mem_mb = 2000
    run:
        import numpy as np
        import tifffile
        import nd2

        with nd2.ND2File(input[0]) as images:
            meta = images.ome_metadata().images[0].dict()
            #debug(meta)

            size_x = meta['pixels']['physical_size_x']# * meta['pixels']['size_x']
            size_y = meta['pixels']['physical_size_y']# * meta['pixels']['size_y']
            xposes = np.array([plane['position_x'] for plane in meta['pixels']['planes'] if plane['the_c'] == 0]) / size_x
            yposes = np.array([plane['position_y'] for plane in meta['pixels']['planes'] if plane['the_c'] == 0]) / size_y
            
            #positions = np.array([-xposes, -yposes]).T
            positions = np.array([yposes, -xposes]).T

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
    input:
        get_nd2filename
    output:
        input_dir + 'well{well}/cycle{cycle}/raw.tif'
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
    input:
        images = expand(input_dir + 'well{well}/cycle{cycle}/raw.tif', cycle=cycles_pt, allow_missing=True),
        positions = expand(input_dir + 'well{well}/cycle{cycle}/positions.csv', cycle=cycles_pt, allow_missing=True),
    output:
        images = expand(input_dir + 'well{well}_subset{size,\d+}/cycle{cycle}/raw.tif', cycle=cycles_pt, allow_missing=True),
        positions = expand(input_dir + 'well{well}_subset{size,\d+}/cycle{cycle}/positions.csv', cycle=cycles_pt, allow_missing=True),
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
                cur_mins = np.round(cur_mins / phenotype_scale)
                cur_maxes = np.round(cur_maxes / phenotype_scale)
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
        low_bound = center - (tile_size - (tile_size // 2))
        high_bound = center + (tile_size // 2)

        for i, poses, path in zip(range(len(all_poses)), all_poses, input.images):
            low, high = low_bound, high_bound
            if any(path.count('cycle' + cycle) for cycle in phenotype_cycles):
                low = low * phenotype_scale
                high = high * phenotype_scale

            images = tifffile.imread(path)
            mask = np.all((low <= poses[:,:2]) & (poses[:,:2] <= high), axis=1)
            debug (path, low, high)
            debug (poses[mask,:2])

            np.savetxt(output.positions[i], poses[mask], delimiter=',', fmt='%d')
            tifffile.imwrite(output.images[i], images[mask])
            del images

