#!/bin/bash
source /net/fowler/vol1/shared/miniconda3/etc/profile.d/conda.sh
conda activate ops 3> /dev/null


mem_arg='$(expr {resources.mem_mb} / {threads})'
out_path='$( (test -f logs/{output[0]}.out && rm logs/{output[0]}.out); mkdir -p $(dirname logs/{output[0]}.out); realpath logs/{output[0]}.out)'
err_path='$( (test -f logs/{output[0]}.err && rm logs/{output[0]}.err); mkdir -p $(dirname logs/{output[0]}.err); realpath logs/{output[0]}.err)'
cuda='$(test "{resources.cuda}" -eq 1 && echo -l cuda=1)'


cluster_cmd="qsub -v PATH,PYTHONNOUSERSITE=1 -terse -l mfree=${mem_arg}M -l h_rt=48:0:0 -l h=fl004 -o $out_path -e $err_path $cuda -pe serial {threads}"

configfile=
if test -f default-config.yaml; then
    configfile='--configfile  default-config.yaml'
fi
if test -f config.yaml; then
    configfile='--configfile  config.yaml'
fi

snakemake \
    --cluster "$cluster_cmd" \
    --cluster-cancel "qdel" \
    -j 128 \
    $* --cores 50 --resources mem_mb=2000000 cuda=2 --set-resource-scopes mem_mb=global threads=global cuda=global \
    --rerun-triggers mtime params input software-env \
    --keep-incomplete --latency-wait 120 \
    --default-resources mem_mb=5000 disk_mb=5000 cuda=0 --use-conda --conda-frontend conda --retries 2


