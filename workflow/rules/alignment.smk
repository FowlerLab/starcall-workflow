import os
import glob

##################################################
##  Aligning tiles and solving for global positions
##################################################

""" Does alignment for one cycle, calculating constraints between neighboring images.
Outputs a composite file, which contains all the calculated data from the cycle
Also outputs plots into the plots directory under plots/scores_well*_cycle*_
which show how well the alignment is going
"""
rule align_cycle:
    input:
        images = stitching_input_dir + '{prefix}/cycle{cycle}/raw.tif',
        rawposes = stitching_input_dir + '{prefix}/cycle{cycle}/positions.csv',
    output:
        composite = stitching_dir + '{prefix}/cycle{cycle}/partial_composite.bin',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb + 5000
    threads: 8
    run:
        import concurrent.futures
        import numpy as np
        import tifffile
        import constitch
        import tifffile


        executor = concurrent.futures.ThreadPoolExecutor(max_workers=threads)
        debug("Running with", threads, "threads")
        debug(input.rawposes, input.images)

        poses = np.loadtxt(input.rawposes, delimiter=',', dtype=int)
        debug ('loading images', wildcards.cycle)

        images = tifffile.memmap(input.images, mode='r')[:,alignment_channel]
        debug(images.shape)

        print_mem('align_poses', resources.mem_mb)
        debug ('stitching cycle', wildcards.cycle)

        composite = constitch.CompositeImage(progress=True, executor=executor)#, precalculate_fft=True)

        input_poses = np.concatenate([poses[:,:2], np.full((len(poses),1), cycles_pt.index(wildcards.cycle), dtype=poses.dtype)], axis=1)
        debug(images.shape, images.ndim)
        if images.ndim == 4:
            composite.add_images(images, input_poses, scale='tile', channel_axis=0)
        else:
            composite.add_images(images, input_poses, scale='tile')
        composite.print_mem_usage()

        overlapping = composite.constraints(touching=True)
        constraints = overlapping.calculate()
        composite.plot_scores('plots/new_scores_{}_cycle{}_step1.png'.format(wildcards.prefix, wildcards.cycle), constraints)

        thresh = composite.calc_score_threshold()
        constraints = constraints.filter(min_score=thresh)
        composite.plot_scores('plots/new_scores_{}_cycle{}_step2.png'.format(wildcards.prefix, wildcards.cycle), constraints)
        stage_model = constraints.fit_model(outliers=True)
        constraints = constraints.filter(stage_model.inliers)
        modeled = overlapping.calculate(stage_model)
        composite.plot_scores('plots/new_scores_{}_cycle{}_step3.png'.format(wildcards.prefix, wildcards.cycle), constraints.merge(modeled))

        if wildcards.cycle in phenotype_cycles:
            composite.boxes.pos1[:,:2] //= phenotype_scale
            composite.boxes.pos2[:,:2] //= phenotype_scale

            for constraint in constraints:
                constraint.dx //= phenotype_scale
                constraint.dy //= phenotype_scale

            for constraint in modeled:
                constraint.dx //= phenotype_scale
                constraint.dy //= phenotype_scale

        debug('done')
        constitch.save(output.composite, composite, constraints, modeled, save_images=False)

""" Combines all the alignment data for each cycle and calculates alignments between the cycles.
outputs the result to one composite file for the whole well, another composite file for each cycle,
and a positions csv file that has the pixel positions for each cycle.
"""
rule align_well:
    input:
        #images = expand(input_dir + '{prefix}/cycle{cycle}.raw.nd2', cycle=cycles_pt, allow_missing=True),
        images = expand(stitching_input_dir + '{prefix}/cycle{cycle}/raw.tif', cycle=cycles_pt, allow_missing=True),
        composites = expand(stitching_dir + '{prefix}/cycle{cycle}/partial_composite.bin', cycle=cycles_pt, allow_missing=True),
        rawposes = expand(stitching_input_dir + '{prefix}/cycle{cycle}/positions.csv', cycle=cycles_pt, allow_missing=True),
    output:
        positions = expand(stitching_dir + '{prefix}/cycle{cycle}/positions.csv', cycle=cycles_pt, allow_missing=True),
        #positions = input_dir + '{prefix}.positions.csv',
        composite = stitching_dir + '{prefix}/composite.bin',
        composites = expand(stitching_dir + '{prefix}/cycle{cycle}/composite.bin', cycle=cycles_pt, allow_missing=True),
    resources:
        mem_mb = lambda wildcards, input: 10000 + input.size_mb / 2.5
    threads: 16
    run:
        import concurrent.futures
        import numpy as np
        import tifffile
        import constitch
        import tifffile

        all_poses = []
        for i,path in enumerate(sorted(input.rawposes)):
            cur_poses = np.loadtxt(path, delimiter=',', dtype=int)
            all_poses.append(cur_poses)
            debug(cur_poses.shape)

        executor = concurrent.futures.ThreadPoolExecutor(max_workers=threads)
        debug("Running with", threads, "threads")
        debug(input.images)

        #"""
        full_composite = constitch.CompositeImage(progress=True, executor=executor)#, precalculate_fft=True)
        all_constraints = constitch.ConstraintSet()
        all_modeled = constitch.ConstraintSet()
        subcomposites = []

        for cycle, rawpath, composite_path in zip(cycles_pt, input.images, input.composites):
            debug ('loading images', cycle)

            debug ('copy')
            images = tifffile.memmap(rawpath, mode='r')[:,alignment_channel].copy()
            #images = tifffile.imread(rawpath)[:,1].copy()
            debug(images.shape)

            composite, constraints, modeled = constitch.load(composite_path)
            composite.images = images

            mean_pos = (composite.boxes.pos1.mean(axis=0) + composite.boxes.pos2.mean(axis=0)) / 2
            if len(subcomposites):
                global_mean_pos = (full_composite.boxes.pos1.mean(axis=0) + full_composite.boxes.pos2.mean(axis=0)) / 2
                debug (mean_pos, global_mean_pos, global_mean_pos - mean_pos)
                offset = np.round(global_mean_pos - mean_pos)[:2].reshape(1,2).astype(int)
                debug (offset)
                composite.boxes.pos1[:,:2] += offset
                composite.boxes.pos2[:,:2] += offset
                debug ((composite.boxes.pos1.mean(axis=0) + composite.boxes.pos2.mean(axis=0)) / 2)

            debug ('composite sizes', len(constraints), len(modeled))
            subcomposite, constraints, modeled = full_composite.merge(composite, constraints, modeled)
            debug ('          sizes', len(constraints), len(modeled))
            debug (' before', len(all_constraints), len(all_modeled))
            subcomposites.append(subcomposite)
            all_constraints.add(constraints)
            all_modeled.add(modeled)
            debug (' after', len(all_constraints), len(all_modeled))

            del images, composite, constraints, modeled

        debug ('adding constraints between cycles')

        max_pairs = 6
        def filter_overlapping(constraint):
            return constraint.overlap > 0 and constraint.box1.pos1[2] != constraint.box2.pos1[2] and abs(constraint.box1.pos1[2] - constraint.box2.pos1[2]) <= max_pairs

        overlapping = full_composite.constraints(filter_overlapping)
        debug ('num pairs', len(overlapping))

        constraints = overlapping.calculate()
        debug (len(constraints), len(all_constraints))
        constraints.add(all_constraints)
        debug (len(constraints))
        debug (constraints[0,1])
        debug ('  done')

        full_composite.plot_scores('plots/new_scores_{}_step1.png'.format(wildcards.prefix), constraints.merge(all_modeled))

        thresh = full_composite.calc_score_threshold()
        debug('thresh', thresh)
        constraints = constraints.filter(min_score=thresh)

        full_composite.plot_scores('plots/new_scores_{}_step2.png'.format(wildcards.prefix), constraints.merge(all_modeled))

        constraints = constraints.filter(max_length=max(full_composite.images[0].shape))

        full_composite.plot_scores('plots/new_scores_{}_step3.png'.format(wildcards.prefix), constraints.merge(all_modeled))

        debug ('Solving constraints')
        full_composite.apply(constraints.merge(all_modeled).solve(constitch.OutlierSolver()))
        #final_poses = full_composite.solve_constraints_old(filter_outliers=True, max_outlier_ratio=0.2)
        #final_poses = full_composite.solve_constraints(filter_outliers=True, max_outlier_ratio=0.2, outlier_threshold=20, scores_plot_path='plots/new_scores_well{}_step3_filters{{}}.png'.format(wildcards.well))
        debug ('  done')

        constitch.save(output.composite, full_composite, constraints, all_modeled, save_images=False)

        full_composite.plot_scores('plots/new_scores_{}_step4.png'.format(wildcards.prefix), constraints.merge(all_modeled))
        full_composite.plot_scores('plots/new_scores_{}_step4_accuracy.png'.format(wildcards.prefix), constraints.merge(all_modeled), score_func='accuracy')

        executor.shutdown()

        for subcomposite, poses, poses_path, composite_path, cycle in zip(subcomposites, all_poses, output.positions, output.composites, cycles_pt):
            if cycle in phenotype_cycles:
                subcomposite.boxes.pos1[:,:2] *= phenotype_scale
                subcomposite.boxes.pos2[:,:2] *= phenotype_scale

            poses = np.concatenate([poses[:,:2], subcomposite.positions], axis=1)
            np.savetxt(poses_path, poses, fmt="%d", delimiter=',')
            constitch.save(composite_path, subcomposite, save_images=False)

