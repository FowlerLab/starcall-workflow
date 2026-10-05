import os
import glob

# cycle names used in the stitching QC plots, indexed by the z plane of the composite
qc_cycle_labels = ['cycle' + cycle for cycle in cycles_pt]

rotation_config = config['stitching'].get('rotation', {})
rotation_order = rotation_config.get('interpolation_order', 1)

def read_cycle_rotations(path):
    """ The rotation (degrees) applied to the tiles of each cycle, from the rotation.csv
    written by detect_rotation, as a dict of cycle name to angle """
    import pandas
    table = pandas.read_csv(path, dtype={'cycle': str})
    return dict(zip(table['cycle'], table['applied_deg'].astype(float)))

##################################################
##  Aligning tiles and solving for global positions
##################################################

rule make_initial_composite:
    """ Creates a constitch.CompositeImage instance with all tiles in the well.
    Each tile is positioned with its tile location from the positions.csv files.
    This provides a base instance with which constraints can be calculated between
    each cycle and between neighboring tiles
    """
    input:
        images = lambda wildcards: [find_input_file(well=wildcards.well_stitching, cycle=cycle) for cycle in cycles_pt],
        #images = expand(input_dir + '{well_stitching}/cycle{cycle}/raw.tif', cycle=cycles_pt, allow_missing=True),
        rawposes = expand(input_dir + '{well_stitching}/cycle{cycle}/positions.csv', cycle=cycles_pt, allow_missing=True),
    output:
        composite = stitching_dir + '{well_stitching}/initial_composite.json',
        plot = qc_dir + '{well_stitching}/initial_composite.png',
    resources:
        mem_mb = 5000
    run:
        import constitch
        import numpy as np

        composite = constitch.CompositeImage()

        for i, cycle in enumerate(cycles_pt):
            subcomposite = composite.layer(i)

            poses = np.loadtxt(input.rawposes[i], delimiter=',', dtype=int)
            poses = poses.reshape(-1, poses.shape[-1])
            poses = poses[:,:2]
            #images = tifffile.memmap(input.images[i], mode='r')[:,0]
            shape, dtype = iminfo(input.images[i])
            if len(shape) == 3:
                shape = (1, *shape)
            images = np.empty(shape[:1] + shape[2:], dtype)
            debug(poses.shape, images.shape)

            subcomposite.add_images(images, poses, scale='tile')
            subcomposite.setimages([None] * len(subcomposite.images))

            # center each cycle
            if config['stitching'].get('center', False):
                center = (subcomposite.boxes.points1.min(axis=0)[:2] + subcomposite.boxes.points2.max(axis=0)[:2]) // 2
                for box in subcomposite.boxes:
                    box.position[:2] -= center

            if cycle in phenotype_cycles:
                for box in subcomposite.boxes:
                    # scaling from phenotype images to base images
                    box.position[:2] *= bases_scale
                    box.position[:2] //= phenotype_scale
                    box.size[:2] *= bases_scale
                    box.size[:2] //= phenotype_scale

            del images

        import starcall.stitching_qc
        starcall.stitching_qc.plot_tile_layout(composite, output.plot, qc_cycle_labels)

        constitch.save(output.composite, composite)

rule detect_rotation:
    """ Measures the rotation of the image content of every cycle on a sample of tiles,
    against the first sequencing cycle (see starcall.rotation), and decides the rotation applied
    to the tiles of each cycle when they are read for alignment and stitching. Cycles rotated
    less than stitching.rotation.threshold degrees from the reference are not rotated.

    rotation.csv has one row per cycle:
        measured_deg: median rotation of the sampled tiles relative to the first sequencing cycle
        relative_deg: rotation relative to the reference (the median of the sequencing cycles,
            or stitching.rotation.reference)
        applied_deg: rotation applied to the tiles of the cycle, 0 if not corrected
        residual_after_deg: rotation relative to the reference measured again after applying
            the correction, should be close to 0
        n_tiles, spread_deg: tiles that could be measured and the interquartile range of their rotation
        scale: scale relative to the reference cycle, for phenotype cycles after rescaling by
            the nominal bases_scale / phenotype_scale
        fit_rms_px: median rms error of the fits, in pixels
        status: rotated, below threshold, not corrected or failed
    rotation_tiles.csv has the measurement of every sampled tile.
    """
    input:
        composite = stitching_dir + '{well_stitching}/initial_composite.json',
        images = lambda wildcards: [find_input_file(well=wildcards.well_stitching, cycle=cycle) for cycle in cycles_pt],
    output:
        table = stitching_dir + '{well_stitching}/rotation.csv',
        tiles = qc_dir + '{well_stitching}/rotation_tiles.csv',
        plot = qc_dir + '{well_stitching}/rotation.png',
    params:
        channel = config['stitching']['channel'],
    threads: config['stitching'].get('threads', 16)
    resources:
        mem_mb = 12000
    run:
        import concurrent.futures
        import pandas
        import constitch
        import starcall.rotation
        import starcall.stitching_qc

        composite = constitch.load(input.composite)
        readers = [starcall.rotation.FrameReader(path, channel_index(params.channel, cycle=cycle))
                for path, cycle in zip(input.images, cycles_pt)]
        cycle_table, tile_table, zero = starcall.rotation.detect_rotations(
                readers, list(cycles_pt), [cycle in phenotype_cycles for cycle in cycles_pt], composite,
                bases_scale / phenotype_scale,
                reference = rotation_config.get('reference', 'median'),
                threshold = rotation_config.get('threshold', 0.05),
                correct = rotation_config.get('correct', True),
                num_tiles = rotation_config.get('num_tiles', 16),
                order = rotation_order,
                executor = concurrent.futures.ThreadPoolExecutor(max_workers=threads))
        for reader in readers:
            reader.close()

        cycle_table, tile_table = pandas.DataFrame(cycle_table), pandas.DataFrame(tile_table)
        debug(cycle_table.to_string())
        cycle_table.to_csv(output.table, index=False)
        tile_table.to_csv(output.tiles, index=False)
        shape, dtype = iminfo(input.images[0])
        starcall.stitching_qc.plot_rotation(cycle_table, tile_table, output.plot,
                rotation_config.get('threshold', 0.05), rotation_config.get('reference', 'median'),
                tile_size=shape[-1], phenotype_cycles=phenotype_cycles)


rule extract_alignment_channel:
    """ Extracts the channel used for alignment from all tiles of a cycle into a
    single channel tif, reading one frame at a time. This means the full multichannel
    image is only read once per cycle, instead of once for every cycle pair in
    calculate_constraints. Tiles are rotated by the rotation detected for the cycle,
    see detect_rotation.
    """
    input:
        images = lambda wildcards: find_input_file(well=wildcards.well_stitching, cycle=wildcards.cycle),
        rotation = stitching_dir + '{well_stitching}/rotation.csv',
    output:
        images = stitching_dir + '{well_stitching}/cycle{cycle}/alignment_channel{channel}.tif',
    params:
        channel = parse_param('channel', config['stitching']['channel']),
    wildcard_constraints:
        channel = '|_channel' + any_channel_regex,
    resources:
        mem_mb = 8000
    run:
        import numpy as np
        import tifffile
        from starcall.rotation import rotate_frame

        chan = channel_index(params.channel, cycle=wildcards.cycle)
        angle = read_cycle_rotations(input.rotation)[wildcards.cycle]
        shape, dtype = iminfo(input.images)
        num_tiles = int(np.prod(shape[:-3])) if len(shape) > 3 else 1
        frame_shape = shape[-2:]

        out = tifffile.memmap(output.images, shape=(num_tiles, *frame_shape), dtype=dtype)

        if input.images.endswith('.nd2'):
            import nd2
            with nd2.ND2File(input.images) as ifile:
                for i in range(num_tiles):
                    out[i] = rotate_frame(ifile.read_frame(i)[chan], angle, rotation_order)
        else:
            images = tifffile.memmap(input.images, mode='r').reshape(-1, *shape[-3:])
            for i in range(num_tiles):
                out[i] = rotate_frame(images[i,chan], angle, rotation_order)
            del images

        out.flush()
        del out


rule calculate_constraints:
    """ Calculates the set of constraints between two cycles, or between adjacent tiles
    in the same cycle if both cycles are the same.
    Each overlapping image between the two cycles is aligned with the phase cross
    correlation algorithm, and the resulting offsets are stored in constraints.json
    In addition a random set of non-overlapping constraints are calculated to estimate
    the score threshold for filtering should be.
    Eg, when loading constraints.json, the following lines could be used:
        composite = constitch.load('{well}/initial_composite.json')
        overlapping, constraints, non_overlapping = constitch.load(
                    '{well_stitching}/initial_composite.json', composite=composite)

    Params:
        channel: the channel used for alignment. Should be an integer index or one
            of the channels listed in the config file.
        subpix: the level of sub pixel precision to calculate, eg 16 would mean
            alignment is done to a 1/16th pixel.
    """
    input:
        composite = stitching_dir + '{well_stitching}/initial_composite.json',
        images1 = stitching_dir + '{well_stitching}/cycle{cycle1}/alignment_channel{channel}.tif',
        images2 = stitching_dir + '{well_stitching}/cycle{cycle2}/alignment_channel{channel}.tif',
        #images1 = input_dir + '{well_stitching}/cycle{cycle1}/raw.tif',
        #images2 = input_dir + '{well_stitching}/cycle{cycle2}/raw.tif',
    output:
        constraints = stitching_dir + '{well_stitching}/cycle{cycle1}/cycle{cycle2}/constraints{channel}{subpix}.json',
        plot = qc_dir + '{well_stitching}/cycle{cycle1}_cycle{cycle2}_scores_calculated{channel}{subpix}.png',
    params:
        channel = parse_param('channel', config['stitching']['channel']),
        subpixel_alignment = parse_param('subpix', config['stitching']['subpixel_alignment']),
    wildcard_constraints:
        channel = '|_channel' + any_channel_regex,
        subpix = '|_subpix\d+',
    resources:
        mem_mb = lambda wildcards, input: size_mb(input) * 1.5 + 5000
    threads: config['stitching'].get('threads', 16)
    run:
        import constitch
        import concurrent.futures
        import numpy as np
        import tifffile

        cycle1, cycle2 = cycles_pt.index(wildcards.cycle1), cycles_pt.index(wildcards.cycle2)

        executor = concurrent.futures.ThreadPoolExecutor(max_workers=max(2, threads))
        composite = constitch.load(input.composite, debug=True, progress=True, executor=executor)
        # images are already the single alignment channel, see extract_alignment_channel.
        # read fully instead of memmapped, random tile access over NFS during alignment was slower
        images = tifffile.imread(input.images1)
        composite.layer(cycle1).setimages(images)

        if cycle1 != cycle2:
            images = tifffile.imread(input.images2)
            composite.layer(cycle2).setimages(images)

            def constraint_filter(const):
                return const.box1.position[2] == cycle1 and const.box2.position[2] == cycle2 and const.overlap_ratio >= 0.1
        else:
            def constraint_filter(const):
                return const.box1.position[2] == cycle1 and const.box2.position[2] == cycle2 and const.touching == True
        debug('images loaded')

        # Only pairs between the two cycles can pass the filter, so only those are checked
        # instead of every pair in the well. Pairs are in the same (i, j) order as
        # composite.pair_func, so the resulting constraints are identical
        layers = composite.boxes.positions[:,2]
        indices1, indices2 = np.flatnonzero(layers == cycle1), np.flatnonzero(layers == cycle2)
        candidate_pairs = [(int(i), int(j)) for i in indices1 for j in indices2 if i < j]
        overlapping = composite.constraints(candidate_pairs).filter(constraint_filter)

        debug ('constraints', len(overlapping), cycle1, cycle2, np.unique(composite.boxes.positions[:,2]))

        calculate_params = {}
        if params.subpixel_alignment != 1:
            calculate_params['aligner'] = constitch.FFTAligner(upscale_factor=params.subpixel_alignment)

        constraints = overlapping.calculate(**calculate_params)

        nonoverlapping = composite.constraints(lambda const:
                const.box1.position[2] == cycle1 and const.box2.position[2] == cycle2
                and const.overlap_x < -3000 and const.overlap_y < -3000, limit=100, random=True)
        erroneous_constraints = nonoverlapping.calculate(**calculate_params)

        import starcall.stitching_qc
        starcall.stitching_qc.plot_pair_constraints(composite, overlapping, constraints, erroneous_constraints, output.plot, qc_cycle_labels)
        constitch.save(output.constraints, overlapping, constraints, erroneous_constraints)

        import resource
        debug('peak rss MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)


rule filter_constraints:
    """ Constraints are filtered using a score threshold, calculated as the 95th percentile
    of the set of non overlapping constraints calcualted in calc_constraints. Any constraints
    with a lower score are removed, and a linear model is fit to the remaining.
    Using the RANSAC algorithm, outliers are removed, and all constraints that were
    removed are replaced with constraints estimated by the linear model
    """
    input:
        composite = stitching_dir + '{well_stitching}/initial_composite.json',
        constraints = stitching_dir + '{well_stitching}/cycle{cycle1}/cycle{cycle2}/constraints{params}.json',
    output:
        constraints = stitching_dir + '{well_stitching}/cycle{cycle1}/cycle{cycle2}/filtered_constraints{params}.json',
        plot = qc_dir + '{well_stitching}/cycle{cycle1}_cycle{cycle2}_scores_filtered{params}.png',
    wildcard_constraints:
        params = params_regex('channel', 'subpix')
    resources:
        mem_mb = lambda wildcards, input: size_mb(input) + 5000
    run:
        import constitch
        import numpy as np

        composite = constitch.load(input.composite)
        overlapping, constraints, erroneous_constraints = constitch.load(input.constraints, composite=composite)

        # kept for the qc plot, which shows which constraints each step removed
        calculated, score_threshold, above_threshold = constraints, None, constraints
        modeled = constitch.ConstraintSet()
        if len(constraints) != 0:
            score_threshold = np.percentile([const.score for const in erroneous_constraints], 95) if len(erroneous_constraints) else 0.5
            constraints = constraints.filter(min_score=score_threshold)
            above_threshold = constraints

            if wildcards.cycle1 == wildcards.cycle2:
                stage_model = constitch.SimpleOffsetModel() if wildcards.cycle1 == wildcards.cycle2 else constitch.GlobalStageModel()
                stage_model = constraints.fit_model(stage_model, outliers=True)
                constraints = stage_model.inliers
                modeled = overlapping.calculate(stage_model)

        import starcall.stitching_qc
        starcall.stitching_qc.plot_filtered_constraints(composite, calculated, erroneous_constraints, score_threshold,
                above_threshold, constraints, modeled, output.plot, qc_cycle_labels)

        constitch.save(output.constraints, constraints, modeled)


def constraints_needed(wildcards):
    cycle_pairs = []
    if wildcards.onlyfirst == '_onlyfirst':
        for i in range(len(cycles_pt)):
            cycle_pairs.append((0, i))
    else:
        for i in range(len(cycles_pt)):
            for j in range(i, min(i + config['stitching'].get('max_cycle_pairs', 16), len(cycles_pt))):
                cycle_pairs.append((i, j))

    paths = []
    for i, j in cycle_pairs:
        paths.append(stitching_dir + '{well_stitching}/' + 'cycle{cycle1}/cycle{cycle2}/filtered_constraints'.format(
            cycle1=cycles_pt[i], cycle2=cycles_pt[j]) + '{params}.json')

    return paths


rule merge_constraints:
    """ All filtered constraints from cycle pairs are combined into composite.json
    Adjusting the stitching.max_cycle_pairs in config.yaml can limit the cycle pairs
    collected, eg if max_cycle_pairs is 5, cycle 0 and cycle 4 would be calculated but
    0 and 5 or 1 and 6 would not.
    """
    input:
        composite = stitching_dir + '{well_stitching}/initial_composite.json',
        constraints = constraints_needed,
    output:
        constraints = stitching_dir + '{well_stitching}/constraints{params}{onlyfirst}.json',
    wildcard_constraints:
        params = params_regex('channel', 'subpix'),
        onlyfirst = '|_onlyfirst',
    run:
        import constitch

        composite = constitch.load(input.composite)

        all_constraints = constitch.ConstraintSet()
        all_modeled = constitch.ConstraintSet()

        for path in input.constraints:
            constraints, modeled = constitch.load(path, composite=composite)
            all_constraints.add(constraints)
            all_modeled.add(modeled)

        constitch.save(output.constraints, all_constraints, all_modeled)


rule solve_constraints:
    """ Finds global positions for each image tile given all constraints.
    Plots are made in the qc dir showing all constraints before and after
    solving.

    Params:
        solver: (mae, mse, spantree, pulp) The type of solver to use.
            mae is default and minimizes mean absolute error. mse minimizes
            mean squared error and spantree constructs a spanning tree.
    """
    input:
        composite = stitching_dir + '{well_stitching}/initial_composite.json',
        constraints = stitching_dir + '{well_stitching}/constraints{params}.json',
    output:
        composite = stitching_dir + '{well_stitching}/composite{params}{solver}.json',
        plot1 = qc_dir + '{well_stitching}/presolve{params}{solver}.png',
        plot2 = qc_dir + '{well_stitching}/solved{params}{solver}.png',
        plot3 = qc_dir + '{well_stitching}/solved_accuracy{params}{solver}.png',
    params:
        solver = parse_param('solver', config['stitching']['solver']),
    wildcard_constraints:
        params = params_regex('channel', 'subpix', 'onlyfirst'),
        solver = '|_solver(mse|mae|spantree|lp|ilp|pulp|rounded)',
    resources:
        mem_mb = lambda wildcards, input: 10000 + size_mb(input) * 500
    threads: lambda wildcards: 8 if wildcards.solver == '_solverpulp' else 1
    run:
        import constitch

        composite = constitch.load(input.composite)

        all_constraints, all_modeled = constitch.load(input.constraints, composite=composite)
        solving_constraints = all_constraints.merge(all_modeled)

        import starcall.stitching_qc
        starcall.stitching_qc.plot_presolve(composite, solving_constraints, output.plot1, qc_cycle_labels)
        initial_positions = composite.boxes.positions.copy()

        if params.solver == 'pulp':
            solution = solving_constraints.solve(solver=params.solver, threads=threads*2)
        else:
            solution = solving_constraints.solve(solver=params.solver)

        composite.setpositions(solution)
        starcall.stitching_qc.plot_solved_positions(composite, initial_positions, output.plot2, qc_cycle_labels)
        starcall.stitching_qc.plot_solved_accuracy(composite, solving_constraints, output.plot3, qc_cycle_labels)

        constitch.save(output.composite, composite)


rule split_composite:
    """ Splits the full well composite up into a composite for a single well
    """
    input:
        composite = stitching_dir + '{well_stitching}/composite{params}.json',
    output:
        composite = stitching_dir + '{well_stitching}/cycle{cycle}/composite{params}.json',
    wildcard_constraints:
        params = params_regex('channel', 'subpix', 'onlyfirst', 'solver', *ashlar_params),
    run:
        import constitch

        cycle = cycles_pt.index(wildcards.cycle)

        composite = constitch.load(input.composite)
        composite = composite.layer(cycle)
        
        if wildcards.cycle in phenotype_cycles:
            for box in composite.boxes:
                # scaling from base images to phenotype
                box.position[:2] *= phenotype_scale
                box.position[:2] //= bases_scale
                box.size[:2] *= phenotype_scale
                box.size[:2] //= bases_scale

        constitch.save(output.composite, composite)



def stitching_qc_filtered_constraints(wildcards):
    """ The filtered_constraints.json of every cycle pair used for composite.json, with default params """
    import types
    return [path.format(well_stitching=wildcards.well_stitching, params='')
            for path in constraints_needed(types.SimpleNamespace(onlyfirst=''))]

rule stitching_qc_plots:
    """ Makes all the stitching QC plots of a well again from files that already exist,
    without recalculating anything: the tile layout, the calculated and filtered constraints
    of every cycle pair, and the plots before and after solving. Useful for wells stitched
    before the plots were changed. To only make the plots, and not rerun any stitching step
    whose inputs look out of date, restrict the run to this rule (--forcerun is needed too,
    otherwise snakemake reports nothing to be done), eg:
        ./run.sh --allowed-rules stitching_qc_plots --forcerun stitching_qc_plots output/qc/well1/stitching_qc
    """
    input:
        initial_composite = stitching_dir + '{well_stitching}/initial_composite.json',
        filtered = stitching_qc_filtered_constraints,
        calculated = lambda wildcards: [path.replace('/filtered_constraints', '/constraints')
                    for path in stitching_qc_filtered_constraints(wildcards)],
        constraints = stitching_dir + '{well_stitching}/constraints.json',
        composite = stitching_dir + '{well_stitching}/composite.json',
    output:
        plots = directory(qc_dir + '{well_stitching}/stitching_qc'),
    resources:
        mem_mb = lambda wildcards, input: 8000 + size_mb(input) * 20
    run:
        import re
        import numpy as np
        import constitch
        import starcall.stitching_qc

        out = output.plots + '/'
        os.makedirs(out, exist_ok=True)

        composite = constitch.load(input.initial_composite)
        starcall.stitching_qc.plot_tile_layout(composite, out + 'initial_composite.png', qc_cycle_labels)

        for calculated_path, filtered_path in zip(input.calculated, input.filtered):
            cycle1, cycle2 = re.search(r'cycle([^/]+)/cycle([^/]+)/filtered_constraints', filtered_path).groups()
            name = out + 'cycle{}_cycle{}_scores_'.format(cycle1, cycle2)
            overlapping, calculated, background = constitch.load(calculated_path, composite=composite)
            starcall.stitching_qc.plot_pair_constraints(composite, overlapping, calculated, background,
                    name + 'calculated.png', qc_cycle_labels)

            # the filtering steps are not saved, so they are repeated as in filter_constraints
            kept, modeled = constitch.load(filtered_path, composite=composite)
            above_threshold, threshold = calculated, None
            if len(calculated):
                threshold = np.percentile([const.score for const in background], 95) if len(background) else 0.5
                above_threshold = calculated.filter(min_score=threshold)
            starcall.stitching_qc.plot_filtered_constraints(composite, calculated, background, threshold,
                    above_threshold, kept, modeled, name + 'filtered.png', qc_cycle_labels)

        all_constraints, all_modeled = constitch.load(input.constraints, composite=composite)
        solving_constraints = all_constraints.merge(all_modeled)
        starcall.stitching_qc.plot_presolve(composite, solving_constraints, out + 'presolve.png', qc_cycle_labels)

        initial_positions = composite.boxes.positions.copy()
        composite.setpositions(constitch.load(input.composite).boxes.positions)
        starcall.stitching_qc.plot_solved_positions(composite, initial_positions, out + 'solved.png', qc_cycle_labels)
        starcall.stitching_qc.plot_solved_accuracy(composite, solving_constraints, out + 'solved_accuracy.png', qc_cycle_labels)
