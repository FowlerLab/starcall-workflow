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
        import fisseq.stitching
        import tifffile


        executor = concurrent.futures.ThreadPoolExecutor(max_workers=threads)
        debug("Running with", threads, "threads")
        debug(input.rawposes, input.images)

        full_composite = fisseq.stitching.CompositeImage(progress=True, executor=executor)

        poses = np.loadtxt(input.rawposes, delimiter=',', dtype=int)
        debug ('loading images', wildcards.cycle)

        images = tifffile.memmap(input.images, mode='r')[:,alignment_channel]
        debug(images.shape)

        print_mem('align_poses', resources.mem_mb)
        debug ('stitching cycle', wildcards.cycle)

        composite = fisseq.stitching.CompositeImage(progress=True, executor=executor)#, precalculate_fft=True)

        input_poses = np.concatenate([poses[:,:2], np.full((len(poses),1), cycles_pt.index(wildcards.cycle), dtype=poses.dtype)], axis=1)
        debug(images.shape, images.ndim)
        if images.ndim == 4:
            composite.add_images(images, input_poses, scale='tile', channel_axis=0)
        else:
            composite.add_images(images, input_poses, scale='tile')
        composite.print_mem_usage()

        composite.calc_constraints()
        composite.plot_scores('plots/scores_{}_cycle{}_step1.png'.format(wildcards.prefix, wildcards.cycle))
        #composite.save('tmp_cycle.bin')

        thresh = composite.calc_score_threshold()
        composite.filter_constraints(thresh)
        composite.plot_scores('plots/scores_{}_cycle{}_step2.png'.format(wildcards.prefix, wildcards.cycle))
        composite.estimate_stage_model(filter_outliers=True)
        composite.model_constraints()
        composite.plot_scores('plots/scores_{}_cycle{}_step3.png'.format(wildcards.prefix, wildcards.cycle))

        if wildcards.cycle in phenotype_cycles:
            composite.boxes.pos1[:,:2] //= phenotype_scale
            composite.boxes.pos2[:,:2] //= phenotype_scale
            for constraint in composite.constraints.values():
                constraint.dx //= phenotype_scale
                constraint.dy //= phenotype_scale

        debug('done')
        composite.save(output.composite, save_images=False)

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
        import fisseq.stitching
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
        full_composite = fisseq.stitching.CompositeImage(progress=True, executor=executor)#, precalculate_fft=True)
        all_indices = []

        for cycle, rawpath, composite_path in zip(cycles_pt, input.images, input.composites):
            debug ('loading images', cycle)

            debug ('copy')
            images = tifffile.memmap(rawpath, mode='r')[:,alignment_channel].copy()
            #images = tifffile.imread(rawpath)[:,1].copy()
            debug(images.shape)

            composite = fisseq.stitching.CompositeImage.load(composite_path)
            composite.images = images

            mean_pos = (composite.boxes.pos1.mean(axis=0) + composite.boxes.pos2.mean(axis=0)) / 2
            if len(all_indices):
                global_mean_pos = (full_composite.boxes.pos1.mean(axis=0) + full_composite.boxes.pos2.mean(axis=0)) / 2
                debug (mean_pos, global_mean_pos, global_mean_pos - mean_pos)
                offset = np.round(global_mean_pos - mean_pos)[:2].reshape(1,2).astype(int)
                debug (offset)
                composite.boxes.pos1[:,:2] += offset
                composite.boxes.pos2[:,:2] += offset
                debug ((composite.boxes.pos1.mean(axis=0) + composite.boxes.pos2.mean(axis=0)) / 2)

            indices = full_composite.merge(composite)
            all_indices.append(indices)
            del images, composite

            #full_composite.align_disconnected_regions()

        debug ('adding constraints between cycles')
        pairs = full_composite.find_unconstrained_pairs(overlap_threshold=[0,0,-100], needs_overlap=True, max_pairs=6)
        debug ('num pairs', len(pairs))

        full_composite.calc_constraints(pairs)
        #full_composite.set_aligner(fisseq.stitching.FFTAligner())
        debug ('  done')
        full_composite.save(output.composite + '.bkp.bin')
        #"""
        #full_composite = fisseq.stitching.CompositeImage.load(output.composite + '.old.bin')
        for pair in list(full_composite.constraints):
            if pair[0] == pair[1]:
                del full_composite.constraints[pair]
        #composite.images = full_composite.images
        #full_composite = composite
        #full_composite.save('tmp1.bin')#, save_images=False)
        
        full_composite.plot_scores('plots/scores_{}_step1.png'.format(wildcards.prefix))

        thresh = full_composite.calc_score_threshold()
        debug('thresh', thresh)
        full_composite.filter_constraints(thresh)
        full_composite.save('tmp2.bin', save_images=False)

        full_composite.plot_scores('plots/scores_{}_step2.png'.format(wildcards.prefix))

        full_composite.filter_outliers(pairs=pairs)

        full_composite.plot_scores('plots/scores_{}_step3.png'.format(wildcards.prefix))

        debug ('Solving constraints')
        full_composite.save('tmp3.bin', save_images=False)
        final_poses = full_composite.solve_constraints_old(filter_outliers=True, max_outlier_ratio=0.2)
        #final_poses = full_composite.solve_constraints(filter_outliers=True, max_outlier_ratio=0.2, outlier_threshold=20, scores_plot_path='plots/scores_well{}_step3_filters{{}}.png'.format(wildcards.well))
        debug ('  done')

        full_composite.save(output.composite, save_images=False)

        full_composite.plot_scores('plots/scores_{}_step4.png'.format(wildcards.prefix))

        executor.shutdown()

        #"""
        #full_composite = fisseq.stitching.CompositeImage.load(input_dir + 'well{}.composite.bin'.format(wildcards.prefix))
        #final_poses = full_composite.boxes.pos1[:,:2]

        for indices, poses, path, composite_path, cycle in zip(all_indices, all_poses, output.positions, output.composites, cycles_pt):
            subcomposite = full_composite.subcomposite(indices)
            if cycle in phenotype_cycles:
                final_poses[indices] *= phenotype_scale
                subcomposite.boxes.pos1[:,:2] *= phenotype_scale
                subcomposite.boxes.pos2[:,:2] *= phenotype_scale

            poses = np.concatenate([poses[:,:2], final_poses[indices]], axis=1)
            np.savetxt(path, poses, fmt="%d", delimiter=',')
            subcomposite.save(composite_path, save_images=False)

