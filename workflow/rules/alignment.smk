import os
import glob

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
        images = expand(input_dir + '{well_stitching}/cycle{cycle}/raw.tif', cycle=cycles_pt, allow_missing=True),
        rawposes = expand(input_dir + '{well_stitching}/cycle{cycle}/positions.csv', cycle=cycles_pt, allow_missing=True),
    output:
        composite = stitching_dir + '{well_stitching}/initial_composite.json',
        plot = qc_dir + '{well_stitching}/initial_composite.png',
    resources:
        mem_mb = 5000
    run:
        import tifffile
        import constitch
        import numpy as np

        composite = constitch.CompositeImage()

        for i, cycle in enumerate(cycles_pt):
            subcomposite = composite.layer(i)

            poses = np.loadtxt(input.rawposes[i], delimiter=',', dtype=int)
            poses = poses[:,:2]
            images = tifffile.memmap(input.images[i], mode='r')[:,0]
            debug(poses.shape, images.shape)

            subcomposite.add_images(images, poses, scale='tile')
            subcomposite.setimages([None] * len(subcomposite.images))

            if cycle in phenotype_cycles:
                for box in subcomposite.boxes:
                    # scaling from phenotype images to base images
                    box.position[:2] *= bases_scale
                    box.position[:2] //= phenotype_scale
                    box.size[:2] *= bases_scale
                    box.size[:2] //= phenotype_scale

            del images

        composite.plot_scores(output.plot, axis_size=24)

        constitch.save(output.composite, composite)

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
        images1 = input_dir + '{well_stitching}/cycle{cycle1}/raw.tif',
        images2 = input_dir + '{well_stitching}/cycle{cycle2}/raw.tif',
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
        mem_mb = lambda wildcards, input: input.size_mb + 5000
    threads: 1
    run:
        import tifffile
        import constitch
        import concurrent.futures
        import numpy as np

        alignment_channel = params.channel

        cycle1, cycle2 = cycles_pt.index(wildcards.cycle1), cycles_pt.index(wildcards.cycle2)

        executor = concurrent.futures.ThreadPoolExecutor(max_workers=max(2, threads))
        composite = constitch.load(input.composite, debug=True, progress=True, executor=executor)
        images = tifffile.memmap(input.images1, mode='r')
        images = images[:,channel_index(alignment_channel,cycle=cycles_pt[cycle1])]

        composite.layer(cycle1).setimages(images)

        if cycle1 != cycle2:
            images = tifffile.memmap(input.images2, mode='r')
            images = images[:,channel_index(alignment_channel,cycle=cycles_pt[cycle2])]
            composite.layer(cycle2).setimages(images)

            def constraint_filter(const):
                return const.box1.position[2] == cycle1 and const.box2.position[2] == cycle2 and const.overlap_ratio >= 0.1
        else:
            def constraint_filter(const):
                return const.box1.position[2] == cycle1 and const.box2.position[2] == cycle2 and const.touching == True

        if cycle1 == cycle2:
            overlapping = composite.constraints(lambda const:
                    const.box1.position[2] == cycle1 and const.box2.position[2] == cycle2 and const.touching == True)
        else:
            overlapping = composite.constraints(lambda const:
                    const.box1.position[2] == cycle1 and const.box2.position[2] == cycle2 and const.overlap_ratio >= 0.1)

        debug ('constraints', len(overlapping), cycle1, cycle2, np.unique(composite.boxes.positions[:,2]))

        calculate_params = {}
        if params.subpixel_alignment != 1:
            calculate_params['aligner'] = constitch.FFTAligner(upscale_factor=params.subpixel_alignment)

        constraints = overlapping.calculate(**calculate_params)

        nonoverlapping = composite.constraints(lambda const:
                const.box1.position[2] == cycle1 and const.box2.position[2] == cycle2
                and const.overlap_x < -3000 and const.overlap_y < -3000, limit=100, random=True)
        erroneous_constraints = nonoverlapping.calculate(**calculate_params)

        composite.plot_scores(output.plot, constraints)
        constitch.save(output.constraints, overlapping, constraints, erroneous_constraints)


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
        mem_mb = lambda wildcards, input: input.size_mb + 5000
    run:
        import constitch
        import numpy as np

        composite = constitch.load(input.composite)
        overlapping, constraints, erroneous_constraints = constitch.load(input.constraints, composite=composite)

        score_threshold = np.percentile([const.score for const in erroneous_constraints], 95) if len(erroneous_constraints) else 0.5
        constraints = constraints.filter(min_score=score_threshold)

        modeled = constitch.ConstraintSet()
        if wildcards.cycle1 == wildcards.cycle2:
            stage_model = constitch.SimpleOffsetModel() if wildcards.cycle1 == wildcards.cycle2 else constitch.GlobalStageModel()
            stage_model = constraints.fit_model(stage_model, outliers=True)
            constraints = stage_model.inliers
            modeled = overlapping.calculate(stage_model)

        composite.plot_scores(output.plot, constraints)

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
        mem_mb = lambda wildcards, input: input.size_mb * 10000 + 25000
    threads: lambda wildcards: 8 if wildcards.solver == '_solverpulp' else 1
    run:
        import constitch

        composite = constitch.load(input.composite)

        all_constraints, all_modeled = constitch.load(input.constraints, composite=composite)
        solving_constraints = all_constraints.merge(all_modeled)

        composite.plot_scores(output.plot1, solving_constraints)

        if params.solver == 'pulp':
            solution = solving_constraints.solve(solver=params.solver, threads=threads*2)
        else:
            solution = solving_constraints.solve(solver=params.solver)

        composite.setpositions(solution)
        composite.plot_scores(output.plot2, solving_constraints)
        composite.plot_scores(output.plot3, solving_constraints, score_func='accuracy')

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



