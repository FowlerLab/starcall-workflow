import os
import sys
import glob


##################################################
## spliting and merging grid
##################################################

def get_composite(wildcards):
    prefix = wildcards.prefix.split('_seqgrid')[0].split('_phenogrid')[0]
    return stitching_output_dir + prefix + '_seqgrid{grid_size}/grid_composite.json'

rule split_grid_table:
    input:
        table = phenotyping_input_dir + '{prefix}/cells.csv',
        composite = get_composite,
    output:
        table = phenotyping_input_dir + '{prefix}_phenogrid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/cells.csv'
    resources:
        mem_mb = lambda wildcards, input: 5000 + input.size_mb * 2
    run:
        import pandas
        import numpy as np
        import fisseq.stitching

        table = pandas.read_csv(input.table, index_col=0)
        composite = fisseq.stitching.CompositeImage.load(input.composite)
        grid_size, x, y = int(wildcards.grid_size), int(wildcards.x), int(wildcards.y)

        box = composite.boxes[x*grid_size+y]
        box.pos1 *= phenotype_scale
        box.pos2 *= phenotype_scale

        contained = []
        for i, cell in table.iterrows():
            cellbox = fisseq.stitching.BBox([cell.bbox_x1, cell.bbox_y1], [cell.bbox_x2, cell.bbox_y2])
            is_contained = box.contains(cellbox)
            if i == 15336: debug(is_contained, cellbox, box)
            for j in range(x*grid_size+y):
                is_contained = is_contained and not composite.boxes[j].contains(cellbox)
            if i == 15336: debug(is_contained, cellbox, box)
            contained.append(is_contained)

        table = table[contained]
        table['bbox_x1'] -= box.pos1[0]
        table['bbox_x2'] -= box.pos1[0]
        table['bbox_y1'] -= box.pos1[1]
        table['bbox_y2'] -= box.pos1[1]
        table['xpos'] -= box.pos1[0]
        table['ypos'] -= box.pos1[1]
        table.to_csv(output.table)

rule split_grid_segmentation:
    input:
        image = phenotyping_input_dir + '{prefix}/{segmentation_type}.tif',
        composite = get_composite,
        table = phenotyping_input_dir + '{prefix}_phenogrid{grid_size}/tile{x}x{y}y/cells.csv',
    output:
        image = phenotyping_input_dir + '{prefix}_phenogrid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y/{segmentation_type}.tif',
    wildcard_constraints:
        segmentation_type = '(cells|nuclei|line|puncta|lines)_mask',
    run:
        import numpy as np
        import tifffile
        import pandas
        import fisseq.stitching
        import skimage.measure

        image = tifffile.memmap(input.image, mode='r')
        table = pandas.read_csv(input.table, index_col=0)

        composite = fisseq.stitching.CompositeImage.load(input.composite)
        grid_size, x, y = int(wildcards.grid_size), int(wildcards.x), int(wildcards.y)

        box = composite.boxes[x*grid_size+y]
        box.pos1 *= phenotype_scale
        box.pos2 *= phenotype_scale

        debug (box)

        section = image[...,box.pos1[0]:box.pos2[0],box.pos1[1]:box.pos2[1]]
        props = skimage.measure.regionprops(section)
        props = {prop.label: prop for prop in props}

        dtype = [dtype for dtype in [np.uint8, np.uint16, np.uint32] if np.iinfo(dtype).max > len(props)][0]
        newimage = np.zeros(section.shape, dtype)

        if wildcards.segmentation_type.count('cells') > 0 or wildcards.segmentation_type.count('nuclei') > 0:
            for newindex, (index, cell) in enumerate(table.iterrows()):
                x1, y1, x2, y2 = int(cell.bbox_x1), int(cell.bbox_y1), int(cell.bbox_x2), int(cell.bbox_y2)
                debug (x1, x2, y1, y2)
                #x1, y1, x2, y2 = cell.bbox_x1 * phenotype_scale, cell.bbox_y1 * phenotype_scale, cell.bbox_x2 * phenotype_scale, cell.bbox_y2 * phenotype_scale
                mask = section[x1:x2,y1:y2] == index
                debug (x1, x2, y1, y2, mask.max(), mask.sum(), index, np.sum(section == index))
                debug (' ', props[index].bbox)
                newimage[x1:x2,y1:y2][mask] = newindex + 1
        else:
            newimage, mapping, reverse_mapping = skimage.segmentation.relabel_sequential(section)
            newimage = newimage.astype(dtype)

        tifffile.imwrite(output.image, newimage)


def get_grid_filenames(wildcards):
    grid_size = int(wildcards.grid_size)
    numbers = ['{:02}'.format(i) for i in range(grid_size)]
    return expand(sequencing_dir + '{possible_output_dir}{prefix}_seqgrid{grid_size}/tile{x}x{y}y/{type}', x=numbers, y=numbers, allow_missing=True)

def get_composite_filename(wildcards):
    #return wildcards.prefix.replace(processing_dir, input_dir) + '_grid{grid_size}/composite.json'
    return wildcards.any_output_dir.replace(output_dir, input_dir) + '{prefix}_grid{grid_size}/grid_composite.json'

rule merge_grid:
    input:
        images = get_grid_filenames,
        composite = stitching_output_dir + '{prefix}_seqgrid{grid_size}/grid_composite.json',
    output:
        image = sequencing_dir + '{possible_output_dir}{prefix}_seqgrid{grid_size,\d+}/{type}',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 10 + 25000
    wildcard_constraints:
        #type = '((?!/tile).)*(csv|tif)',
        type = '[^/]+',
        possible_output_dir = '(' + output_dir + ')|',
        #any_processing_dir = sequencing_dir + '|' + phenotyping_dir
    run:
        import numpy as np
        import tifffile
        import fisseq.stitching
        import pandas

        composite = fisseq.stitching.CompositeImage.load(input.composite)

        if wildcards.type == 'bases.csv':
            filtered = []
            for i,path in enumerate(input.images):
                bases = np.loadtxt(path, delimiter=',')
                bases[:,:2] += composite.boxes[i].pos1[:2].reshape(1,2)
                filtered.append(bases)

            np.savetxt(output.image, np.concatenate(filtered, axis=0), delimiter=',', fmt='%f')

        elif wildcards.type.endswith('.tif'):

            for i,path in enumerate(input.images):
                composite.images[i] = tifffile.imread(path)

            merger = fisseq.stitching.LastMerger()
            if wildcards.type.endswith('_mask_downscaled.tif'):
                merger = fisseq.stitching.MaskMerger()
            if wildcards.type.endswith('_mask.tif'):
                composite.boxes.pos1 *= phenotype_scale
                composite.boxes.pos2 *= phenotype_scale
                merger = fisseq.stitching.MaskMerger()

            full_image = composite.stitch_images(merger=merger)
            del composite
            tifffile.imwrite(output.image, full_image)


rule link_seq_grid_input:
    input:
        stitching_output_dir + '{prefix}_seqgrid{grid_size}/tile{x}x{y}y{possible_cycle}/{corrected}_pt.tif',
    output:
        phenotyping_input_dir + '{prefix}{possible_seqgrid}_phenogrid{grid_size,\d+}/tile{x,\d+}x{y,\d+}y{possible_cycle}/{corrected,raw|corrected}_pt.tif',
    wildcard_constraints:
        prefix = '((?!_seqgrid).)*',
        possible_seqgrid = '(_seqgrid\d+)|',
        possible_cycle = '(/cycle(' + '|'.join(cycles_pt) + '))|',
    localrule: True
    run:
        os.system('ln -s "{}" "{}"'.format(os.path.relpath(input[0], os.path.dirname(output[0])), output[0]))

rule link_grid:
    input:
        stitching_output_dir + '{prefix}/{type}'
    output:
        phenotyping_input_dir + '{prefix}_seqgrid{grid_size,\d+}/{type}'
    resources:
        mem_mb = 1000
    wildcard_constraints:
        type = 'raw_pt.tif|corrected_pt.tif'
    localrule: True
    run:
        os.system('ln -s "{}" "{}"'.format(os.path.relpath(input[0], os.path.dirname(output[0])), output[0]))

ruleorder: link_grid > link_input_phenotyping


