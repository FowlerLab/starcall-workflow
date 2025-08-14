import os
import glob

alignment_channel = config.get('alignment_channel', 0)
max_constraint_pairs = config.get('max_constraint_pairs', 9999)

##################################################
##  Aligning tiles and solving for global positions
##################################################

rule make_initial_composite:
    input:
        images = expand(stitching_input_dir + '{prefix}/cycle{cycle}/raw.tif', cycle=cycles_pt, allow_missing=True),
        rawposes = expand(stitching_input_dir + '{prefix}/cycle{cycle}/positions.csv', cycle=cycles_pt, allow_missing=True),
    output:
        composite = stitching_dir + '{prefix}/initial_composite.json',
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
            images = tifffile.memmap(input.images[i], mode='r')[:,alignment_channel]
            debug(poses.shape, images.shape)

            subcomposite.add_images(images, poses, scale='tile')
            subcomposite.setimages([None] * len(subcomposite.images))

            if cycle in phenotype_cycles:
                for box in subcomposite.boxes:
                    box.position[:2] //= phenotype_scale
                    box.size[:2] //= phenotype_scale

            del images

        constitch.save(output.composite, composite)

rule calculate_constraints:
    input:
        composite = stitching_dir + '{prefix}/initial_composite.json',
        images1 = stitching_input_dir + '{prefix}/cycle{cycle1}/raw.tif',
        images2 = stitching_input_dir + '{prefix}/cycle{cycle2}/raw.tif',
    output:
        constraints = stitching_dir + '{prefix}/cycle{cycle1}/cycle{cycle2}/constraints.json',
        plot = qc_dir + '{prefix}/cycle{cycle1}_cycle{cycle2}_scores_calculated.png',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb + 5000
    threads: 1
    run:
        import tifffile
        import constitch
        import concurrent.futures
        import numpy as np

        cycle1, cycle2 = cycles_pt.index(wildcards.cycle1), cycles_pt.index(wildcards.cycle2)

        executor = concurrent.futures.ThreadPoolExecutor(max_workers=max(2, threads))
        composite = constitch.load(input.composite, debug=True, progress=True, executor=executor)
        images = tifffile.memmap(input.images1, mode='r')[:,alignment_channel]
        composite.layer(cycle1).setimages(images)

        if cycle1 != cycle2:
            images = tifffile.memmap(input.images2, mode='r')[:,alignment_channel]
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
        #constraints = overlapping.calculate()
        constraints = overlapping.calculate(constitch.FFTAligner(upscale_factor=16))
        composite.plot_scores('plots/tmp_scores.png', constraints)

        nonoverlapping = composite.constraints(lambda const:
                const.box1.position[2] == cycle1 and const.box2.position[2] == cycle2
                and const.overlap_x < -3000 and const.overlap_y < -3000, limit=100, random=True)
        erroneous_constraints = nonoverlapping.calculate()

        composite.plot_scores(output.plot, constraints)

        constitch.save(output.constraints, overlapping, constraints, erroneous_constraints)


rule filter_constraints:
    input:
        composite = stitching_dir + '{prefix}/initial_composite.json',
        constraints = stitching_dir + '{prefix}/cycle{cycle1}/cycle{cycle2}/constraints.json',
    output:
        constraints = stitching_dir + '{prefix}/cycle{cycle1}/cycle{cycle2}/filtered_constraints.json',
        plot = qc_dir + '{prefix}/cycle{cycle1}_cycle{cycle2}_scores_filtered.png',
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
    paths = []
    for i in range(len(cycles_pt)):
        for j in range(i, min(i + max_constraint_pairs, len(cycles_pt))):
            paths.append(stitching_dir + '{prefix}/' + 'cycle{cycle1}/cycle{cycle2}/filtered_constraints.json'.format(
                cycle1=cycles_pt[i], cycle2=cycles_pt[j]))
    return paths

rule merge_constraints:
    input:
        composite = stitching_dir + '{prefix}/initial_composite.json',
        constraints = constraints_needed,
    output:
        constraints = stitching_dir + '{prefix}/constraints.json',
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
    input:
        composite = stitching_dir + '{prefix}/initial_composite.json',
        constraints = stitching_dir + '{prefix}/constraints.json',
    output:
        composite = stitching_dir + '{prefix,[^/]*}/composite.json',
        plot1 = qc_dir + '{prefix}/presolve.png',
        plot2 = qc_dir + '{prefix}/solved.png',
        plot3 = qc_dir + '{prefix}/solved_accuracy.png',
    resources:
        mem_mb = lambda wildcards, input: input.size_mb * 1000 + 25000
    run:
        import constitch

        composite = constitch.load(input.composite)

        all_constraints, all_modeled = constitch.load(input.constraints, composite=composite)
        solving_constraints = all_constraints.merge(all_modeled)

        composite.plot_scores(output.plot1, solving_constraints)

        solution = solving_constraints.solve(solver='mae')

        composite.setpositions(solution)
        composite.plot_scores(output.plot2, solving_constraints)
        composite.plot_scores(output.plot3, solving_constraints, score_func='accuracy')

        constitch.save(output.composite, composite)


rule split_composite:
    input:
        composite = stitching_dir + '{prefix}/composite.json',
    output:
        composite = stitching_dir + '{prefix}/cycle{cycle}/composite.json',
    run:
        import constitch

        cycle = cycles_pt.index(wildcards.cycle)

        composite = constitch.load(input.composite)
        composite = composite.layer(cycle)
        
        if wildcards.cycle in phenotype_cycles:
            for box in composite.boxes:
                box.position[:2] *= phenotype_scale
                box.size[:2] *= phenotype_scale

        constitch.save(output.composite, composite)



