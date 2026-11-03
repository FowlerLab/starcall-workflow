import os
import glob

# cycle names used in the stitching QC plots, indexed by the z plane of the composite
qc_cycle_labels = ['cycle' + cycle for cycle in cycles_pt]

rotation_config = config['stitching'].get('rotation', {})
rotation_order = rotation_config.get('interpolation_order', 1)
# 'resample' rotates and rescales the pixels of phenotype tiles in every stitched image, 'coordinates'
# only for alignment, and keeps the correction as a transform, see starcall.phenotype_frame
phenotype_correction = rotation_config.get('phenotype_correction', 'resample')
assert phenotype_correction in ('resample', 'coordinates'), 'stitching.rotation.phenotype_correction should be resample or coordinates'
phenotype_coordinates = phenotype_correction == 'coordinates'
phenotype_frame_config = config['stitching'].get('phenotype_frame', {})

def read_cycle_rotations(path):
    """ The rotation (degrees) applied to the tiles of each cycle, from the rotation.csv
    written by detect_rotation, as a dict of cycle name to angle """
    import pandas
    table = pandas.read_csv(path, dtype={'cycle': str})
    return dict(zip(table['cycle'], table['applied_deg'].astype(float)))

def read_cycle_scales(path):
    """ The scale applied to the content of the tiles of each cycle, from rotation.csv,
    as a dict of cycle name to scale (1 for every cycle in tables from before scale correction) """
    import pandas
    table = pandas.read_csv(path, dtype={'cycle': str})
    if 'applied_scale' not in table:
        return {cycle: 1.0 for cycle in table['cycle']}
    return dict(zip(table['cycle'], table['applied_scale'].astype(float)))

def stitch_corrections(path, stitched_cycles):
    """ (angles, scales) the tiles of each of stitched_cycles are resampled by when stitched, from rotation.csv.
    With phenotype_correction: coordinates, phenotype cycles are only resampled by the difference between
    their correction and the first phenotype cycle's (none for that one), see starcall.phenotype_frame """
    angles, scales = read_cycle_rotations(path), read_cycle_scales(path)
    if phenotype_coordinates:
        reference = str(phenotype_cycles[0])
        for cycle in phenotype_cycles:
            if str(cycle) != reference:
                angles[cycle] -= angles[reference]
                scales[cycle] /= scales[reference]
        angles[reference], scales[reference] = 0.0, 1.0
    return [angles[cycle] for cycle in stitched_cycles], [scales[cycle] for cycle in stitched_cycles]

##################################################
##  Aligning tiles and solving for global positions
##################################################

rule make_initial_composite:
    """ Creates a constitch.CompositeImage instance with all tiles in the well.
    Each tile is positioned with its tile location from the positions.csv files.
    This provides a base instance with which constraints can be calculated between
    each cycle and between neighboring tiles

    If stitching.phenotype_stage_positions is set, the phenotype cycles are instead placed
    from their stage positions relative to the first sequencing cycle (see starcall.stage),
    this should be used only for a phenotype cycle imaged as a single tile that doesn't cover the whole well.
    """
    input:
        images = lambda wildcards: [find_input_file(well=wildcards.well_stitching, cycle=cycle) for cycle in cycles_pt],
        #images = expand(input_dir + '{well_stitching}/cycle{cycle}/raw.tif', cycle=cycles_pt, allow_missing=True),
        rawposes = expand(input_dir + '{well_stitching}/cycle{cycle}/positions.csv', cycle=cycles_pt, allow_missing=True),
        stage = lambda wildcards: expand(input_dir + '{well}/cycle{cycle}/stage.csv', well=wildcards.well_stitching, cycle=cycles_pt)
                if config['stitching'].get('phenotype_stage_positions', False) else [],
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

        if input.stage:
            import starcall.stage
            stages = [np.loadtxt(path, delimiter=',', ndmin=2) for path in input.stage]
            reference = [i for i, cycle in enumerate(cycles_pt) if cycle not in phenotype_cycles][0]
            moved = [i for i, cycle in enumerate(cycles_pt) if cycle in phenotype_cycles]
            rms = starcall.stage.place_from_stage(composite, stages, reference, moved)
            debug('stage map fit to cycle', cycles_pt[reference], 'rms px', rms)
            for i in moved:
                debug('cycle', cycles_pt[i], 'placed from stage positions, center',
                        starcall.stage.layer_centers(composite, i).mean(axis=0))

        import starcall.stitching_qc
        starcall.stitching_qc.plot_tile_layout(composite, output.plot, qc_cycle_labels)

        constitch.save(output.composite, composite)

rule detect_rotation:
    """ Measures the rotation of the image content of every cycle on a sample of tiles,
    against the first sequencing cycle (see starcall.rotation), and decides the rotation applied
    to the tiles of each cycle when they are read for alignment and stitching. Cycles rotated
    less than stitching.rotation.threshold degrees from the reference are not rotated.
    A phenotype cycle imaged as a single tile is measured on the tiles of the first sequencing
    cycle inside it (at least stitching.rotation.min_phenotype_overlap of their area), instead of
    the tiles sampled over the well.

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
        scale_spread: interquartile range of the scale of the tiles
        applied_scale: scale applied to the tiles of the cycle (stitching.rotation.correct_scale,
            phenotype cycles only), 1 if not corrected
        scale_after: scale measured again after applying the correction, should be close to 1
        fit_rms_px: median rms error of the fits, in pixels
        status: rotated, below threshold, not corrected or failed, and rescaled
    rotation_tiles.csv has the measurement of every sampled tile. rotation.png has a summary of
    the correction of each cycle and how it is applied (stitching.rotation.phenotype_correction).
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
                executor = concurrent.futures.ThreadPoolExecutor(max_workers=threads),
                min_phenotype_overlap = rotation_config.get('min_phenotype_overlap', 0.5),
                correct_scale = rotation_config.get('correct_scale', False),
                scale_threshold = rotation_config.get('scale_threshold', 0.0005))
        for reader in readers:
            reader.close()

        cycle_table, tile_table = pandas.DataFrame(cycle_table), pandas.DataFrame(tile_table)
        debug(cycle_table.to_string())
        cycle_table.to_csv(output.table, index=False)
        tile_table.to_csv(output.tiles, index=False)
        shape, dtype = iminfo(input.images[0])

        # the linear part of the transform of the phenotype frame, the translation is only known after solving
        frame, pt_shape = None, None
        if phenotype_cycles:
            import starcall.phenotype_frame
            frame = starcall.phenotype_frame.make_frame(cycle_table, phenotype_cycles, phenotype_scale / bases_scale, (0, 0), (0, 0))
            pt_shape = iminfo(input.images[cycles_pt.index(phenotype_cycles[0])])[0][-2:]
        summary = starcall.stitching_qc.rotation_summary(cycle_table, phenotype_cycles, phenotype_scale / bases_scale,
                phenotype_correction, tile_size=pt_shape, frame=frame)
        debug('\n'.join(summary))
        starcall.stitching_qc.plot_rotation(cycle_table, tile_table, output.plot,
                rotation_config.get('threshold', 0.05), rotation_config.get('reference', 'median'),
                tile_size=shape[-1], phenotype_cycles=phenotype_cycles, summary=summary)


rule extract_alignment_channel:
    """ Extracts the channel used for alignment from all tiles of a cycle into a
    single channel tif, reading one frame at a time. This means the full multichannel
    image is only read once per cycle, instead of once for every cycle pair in
    calculate_constraints. Tiles are rotated and rescaled by the correction detected for the cycle,
    see detect_rotation. This is done for phenotype cycles too when they are corrected in coordinates
    (stitching.rotation.phenotype_correction), so that the solved positions are those of corrected tiles.
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
        scale = read_cycle_scales(input.rotation)[wildcards.cycle]
        shape, dtype = iminfo(input.images)
        num_tiles = int(np.prod(shape[:-3])) if len(shape) > 3 else 1
        frame_shape = shape[-2:]

        out = tifffile.memmap(output.images, shape=(num_tiles, *frame_shape), dtype=dtype)

        if input.images.endswith('.nd2'):
            import nd2
            with nd2.ND2File(input.images) as ifile:
                for i in range(num_tiles):
                    out[i] = rotate_frame(ifile.read_frame(i)[chan], angle, rotation_order, scale)
        else:
            images = tifffile.memmap(input.images, mode='r').reshape(-1, *shape[-3:])
            for i in range(num_tiles):
                out[i] = rotate_frame(images[i,chan], angle, rotation_order, scale)
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

    If stitching.phenotype_cropped_aligner is set, pairs between a sequencing and a phenotype
    cycle are instead aligned with starcall.cropped_alignment.CroppedFFTAligner, which aligns each
    sequencing tile against a window of the phenotype tile around its expected position, within
    stitching.phenotype_cropped_margin px. No non overlapping constraints are calculated for these.

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
        mem_mb = lambda wildcards, input, attempt: (size_mb(input) * 1.5 + 5000) * attempt
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

        phenotype_pair = (wildcards.cycle1 in phenotype_cycles) != (wildcards.cycle2 in phenotype_cycles)
        if phenotype_pair and config['stitching'].get('phenotype_cropped_aligner', False):
            import starcall.cropped_alignment
            aligner = starcall.cropped_alignment.CroppedFFTAligner(
                    margin=config['stitching'].get('phenotype_cropped_margin', 600),
                    upscale_factor=params.subpixel_alignment)
            constraints = overlapping.calculate(aligner=aligner)
            # the cropped aligner only searches near the expected offset, so the score of pairs
            # far apart doesn't estimate the score of erroneous constraints
            erroneous_constraints = constitch.ConstraintSet()
        else:
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

    If stitching.phenotype_score_threshold is set, constraints between a sequencing and a
    phenotype cycle are instead filtered with that score threshold, and the remaining are
    checked for consistency: constraints whose correction from the initial positions is more than
    stitching.phenotype_consistency_tolerance px from an affine fit of the others are removed,
    see starcall.cropped_alignment.consistent_constraints
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
        phenotype_pair = (wildcards.cycle1 in phenotype_cycles) != (wildcards.cycle2 in phenotype_cycles)
        phenotype_threshold = config['stitching'].get('phenotype_score_threshold', None)
        if len(constraints) != 0 and phenotype_pair and phenotype_threshold is not None:
            import starcall.cropped_alignment
            score_threshold = phenotype_threshold
            constraints = constraints.filter(min_score=score_threshold)
            above_threshold = constraints
            constraints = starcall.cropped_alignment.consistent_constraints(constraints,
                    tolerance=config['stitching'].get('phenotype_consistency_tolerance', 20))
            debug('phenotype pair: {} calculated, {} above score threshold {}, {} consistent'.format(
                    len(calculated), len(above_threshold), score_threshold, len(constraints)))

        elif len(constraints) != 0:
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


rule calculate_frame_constraints:
    """ Aligns the overlapping tiles of a phenotype cycle with each other at full phenotype resolution,
    on the tiles as they are stitched in the phenotype frame (raw frames for the first phenotype cycle),
    for phenotype_frame. Only used with stitching.rotation.phenotype_correction: coordinates.
    Tiles are placed at their frame positions from the joint solve, see starcall.phenotype_frame.
    Constraints are filtered like the constraints between tiles of a cycle in filter_constraints:
    score threshold from non overlapping pairs, then a SimpleOffsetModel fit whose outliers are removed
    and which estimates the pairs without a constraint.

    frame_constraints.csv has one row per pair: index1, index2 (tile index within the cycle), dx, dy
    (offset in phenotype px from tile index1 to tile index2, row and column), score, and
    kind: measured or modeled. It is empty for a cycle imaged as a single tile.
    """
    input:
        composite = stitching_dir + '{well_stitching}/composite.json',
        rotation = stitching_dir + '{well_stitching}/rotation.csv',
        images = lambda wildcards: find_input_file(well=wildcards.well_stitching, cycle=wildcards.cycle),
    output:
        table = stitching_dir + '{well_stitching}/cycle{cycle}/frame_constraints.csv',
        plot = qc_dir + '{well_stitching}/cycle{cycle}_frame_constraints.png',
    params:
        channel = config['stitching']['channel'],
        subpixel_alignment = config['stitching'].get('subpixel_alignment', 1),
    wildcard_constraints:
        cycle = '|'.join(phenotype_cycles) if phenotype_cycles else 'PT',
    threads: config['stitching'].get('threads', 16)
    resources:
        mem_mb = lambda wildcards, input, attempt: (size_mb(input) * 1.5 + 8000) * attempt
    run:
        import concurrent.futures
        import numpy as np
        import pandas
        import constitch
        import starcall.phenotype_frame
        import starcall.stitching_qc
        from starcall.rotation import FrameReader, rotate_frame

        cycle = cycles_pt.index(wildcards.cycle)
        full = constitch.load(input.composite, constraints=False)
        mins, maxes = full.boxes.points1.min(axis=0)[:2], full.boxes.points2.max(axis=0)[:2]
        frame = starcall.phenotype_frame.make_frame(pandas.read_csv(input.rotation, dtype={'cycle': str}),
                phenotype_cycles, phenotype_scale / bases_scale, mins, maxes)
        layer = full.layer(cycle)
        n_tiles = len(layer.boxes)
        columns = ['index1', 'index2', 'dx', 'dy', 'score', 'kind']

        if n_tiles < 2:
            debug('cycle', wildcards.cycle, 'has', n_tiles, 'tile, no constraints between its tiles')
            pandas.DataFrame(columns=columns).to_csv(output.table, index=False)
            starcall.stitching_qc.plot_tile_layout(full, output.plot, qc_cycle_labels,
                    title='cycle{} is a single tile, no constraints between its tiles'.format(wildcards.cycle))
        else:
            angle, scale = frame.relative(wildcards.cycle)
            reader = FrameReader(input.images, channel_index(params.channel, cycle=wildcards.cycle))
            images = np.array([rotate_frame(reader(i), angle, rotation_order, scale) for i in range(n_tiles)])
            reader.close()
            sizes = np.tile(images.shape[1:3], (n_tiles, 1))
            anchors = frame.anchor_positions(wildcards.cycle, layer.boxes.positions[:,:2], sizes)

            executor = concurrent.futures.ThreadPoolExecutor(max_workers=max(2, threads))
            composite = constitch.CompositeImage(executor=executor, debug=True, progress=True)
            composite.add_images(images, np.rint(anchors).astype(int))

            # as calculate_constraints for a cycle with itself
            overlapping = composite.constraints(lambda const: const.touching == True)
            aligner = constitch.FFTAligner(upscale_factor=params.subpixel_alignment) if params.subpixel_alignment != 1 else 'fft'
            calculated = overlapping.calculate(aligner=aligner)
            # far apart pairs give the score of a wrong alignment, -3000 base px as in calculate_constraints
            far = -3000 * phenotype_scale / bases_scale
            erroneous = composite.constraints(lambda const: const.overlap_x < far and const.overlap_y < far,
                    limit=100, random=True).calculate(aligner=aligner)

            # as filter_constraints
            threshold = np.percentile([const.score for const in erroneous], 95) if len(erroneous) else 0.5
            above_threshold = calculated.filter(min_score=threshold)
            if len(above_threshold):
                stage_model = above_threshold.fit_model(constitch.SimpleOffsetModel(), outliers=True)
                kept, modeled = stage_model.inliers, overlapping.calculate(stage_model)
            else:
                kept, modeled = above_threshold, constitch.ConstraintSet()
            debug('cycle {}: {} overlapping pairs, {} calculated, {} above score threshold {:.3f}, {} kept, {} modeled'.format(
                    wildcards.cycle, len(overlapping), len(calculated), len(above_threshold), threshold, len(kept), len(modeled)))
            starcall.stitching_qc.plot_filtered_constraints(composite, calculated, erroneous, threshold,
                    above_threshold, kept, modeled, output.plot, ['cycle' + wildcards.cycle])

            rows = [dict(index1=c.index1, index2=c.index2, dx=c.dx, dy=c.dy, score=c.score, kind='measured') for c in kept]
            measured = set((c.index1, c.index2) for c in kept)
            rows += [dict(index1=c.index1, index2=c.index2, dx=c.dx, dy=c.dy, score=np.nan, kind='modeled')
                    for c in modeled if (c.index1, c.index2) not in measured]
            pandas.DataFrame(rows, columns=columns).to_csv(output.table, index=False)


rule phenotype_frame:
    """ The phenotype frame of a well, used with stitching.rotation.phenotype_correction: coordinates
    (see starcall.phenotype_frame): the transform between the stitched phenotype images and the sequencing
    images, from the rotation and scale in rotation.csv, and the positions of the tiles of each phenotype
    cycle in it, where their raw frames are placed by translation only.

    The positions follow from the joint solve (composite.json). For phenotype cycles with more than one tile,
    they are refined with the offsets between overlapping tiles measured at full resolution
    (calculate_frame_constraints), softly anchored to the joint solve:
        stitching.phenotype_frame.refine_weight: weight holding each tile to its position from the joint
            solve, relative to a measured pair (1.0); larger values follow the joint solve more closely
        stitching.phenotype_frame.modeled_weight: weight of the pairs estimated by the stage model (0.25)
        stitching.phenotype_frame.outlier_px: pairs disagreeing with the joint solve by more than this are
            dropped, unset to keep all

    phenotype_frame.json has the transform (frame_to_base: base = A @ frame + t, in (row, col) global composite
    coordinates) and the integer positions of the tiles of each phenotype cycle. phenotype_frame.txt and
    phenotype_frame.png describe it.
    """
    input:
        composite = stitching_dir + '{well_stitching}/composite.json',
        rotation = stitching_dir + '{well_stitching}/rotation.csv',
        constraints = expand(stitching_dir + '{well_stitching}/cycle{cycle}/frame_constraints.csv', cycle=phenotype_cycles, allow_missing=True),
        images = lambda wildcards: [find_input_file(well=wildcards.well_stitching, cycle=cycle) for cycle in phenotype_cycles],
    output:
        frame = stitching_dir + '{well_stitching}/phenotype_frame.json',
        text = qc_dir + '{well_stitching}/phenotype_frame.txt',
        plot = qc_dir + '{well_stitching}/phenotype_frame.png',
    resources:
        mem_mb = 8000
    run:
        import numpy as np
        import pandas
        import constitch
        import starcall.phenotype_frame
        import starcall.stitching_qc

        full = constitch.load(input.composite, constraints=False)
        mins, maxes = full.boxes.points1.min(axis=0)[:2], full.boxes.points2.max(axis=0)[:2]
        cycle_table = pandas.read_csv(input.rotation, dtype={'cycle': str})
        frame = starcall.phenotype_frame.make_frame(cycle_table, phenotype_cycles, phenotype_scale / bases_scale, mins, maxes)
        lam = phenotype_frame_config.get('refine_weight', 1.0)
        modeled_weight = phenotype_frame_config.get('modeled_weight', 0.25)
        outlier_px = phenotype_frame_config.get('outlier_px', None)

        anchors, refined_positions, positions, sizes, stats = {}, {}, {}, {}, {}
        for cycle, constraints_path, images_path in zip(phenotype_cycles, input.constraints, input.images):
            layer = full.layer(cycles_pt.index(cycle))
            shape = iminfo(images_path)[0][-2:]
            sizes[cycle] = np.tile(shape, (len(layer.boxes), 1))
            anchors[cycle] = frame.anchor_positions(cycle, layer.boxes.positions[:,:2], sizes[cycle])
            table = pandas.read_csv(constraints_path)
            weights = np.where(table['kind'] == 'measured', 1.0, modeled_weight)
            refined, used = starcall.phenotype_frame.refine_positions(anchors[cycle], table['index1'], table['index2'],
                    table[['dx', 'dy']].to_numpy(float), weights, lam=lam, outlier_px=outlier_px)
            refined_positions[cycle], positions[cycle] = refined, np.rint(refined).astype(int)
            moved = np.linalg.norm(refined - anchors[cycle], axis=1)
            stats[cycle] = dict(tiles=len(refined), pairs=int(len(table)), pairs_used=int(used),
                    moved_median_px=float(np.median(moved)), moved_max_px=float(moved.max()))

        starcall.phenotype_frame.save(output.frame, frame, positions, extra=dict(refinement=stats,
                refine_weight=lam, modeled_weight=modeled_weight, outlier_px=outlier_px))

        shape = iminfo(input.images[0])[0][-2:]
        summary = ['Phenotype frame of {}: raw phenotype tiles are placed by translation, the rotation and scale'.format(
                wildcards.well_stitching), 'are applied in coordinates when phenotype masks are matched to sequencing dots.', '']
        summary += starcall.stitching_qc.rotation_summary(cycle_table, phenotype_cycles, phenotype_scale / bases_scale,
                phenotype_correction)
        summary += [''] + starcall.phenotype_frame.describe(frame, shape) + ['']
        for cycle, stat in stats.items():
            summary.append('cycle{}: {} tiles, {} pairs between tiles ({} used), refinement moved tiles median {:.2f} max {:.2f} phenotype px'.format(
                    cycle, stat['tiles'], stat['pairs'], stat['pairs_used'], stat['moved_median_px'], stat['moved_max_px']))
        debug('\n'.join(summary))
        with open(output.text, 'w') as ofile:
            ofile.write('\n'.join(summary) + '\n')
        starcall.stitching_qc.plot_phenotype_frame(frame, anchors, refined_positions, sizes, output.plot, summary)


rule split_composite:
    """ Splits the full well composite up into a composite for a single well.
    Phenotype cycles are scaled to phenotype pixels. With stitching.rotation.phenotype_correction:
    coordinates, their tiles are instead placed at their positions in the phenotype frame (phenotype_frame).
    """
    input:
        composite = stitching_dir + '{well_stitching}/composite{params}.json',
        frame = lambda wildcards: stitching_dir + '{well_stitching}/phenotype_frame.json'
                if phenotype_coordinates and wildcards.cycle in phenotype_cycles and wildcards.params == '' else [],
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

        if input.frame:
            import starcall.phenotype_frame
            frame, positions = starcall.phenotype_frame.load(input.frame)
            assert len(positions[wildcards.cycle]) == len(composite.boxes)
            for box, position in zip(composite.boxes, positions[wildcards.cycle]):
                box.position[:2] = position

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
