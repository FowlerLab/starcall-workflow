#!/bin/bash

conda activate ops

test -e logs || mkdir logs
test -e logs/latest && rm -r logs/latest
mkdir logs/latest

mem_arg='$(expr {resources.mem_mb} / {threads})'
#out_path=$(realpath ./logs/latest/)
#err_path=$(realpath ./logs/latest/)
#out_path='$(realpath {log})'
#err_path='$(realpath {log})'
out_path='$( (test -f logs/{output[0]}.out && rm logs/{output[0]}.out); mkdir -p $(dirname logs/{output[0]}.out); realpath logs/{output[0]}.out)'
err_path='$( (test -f logs/{output[0]}.err && rm logs/{output[0]}.err); mkdir -p $(dirname logs/{output[0]}.err); realpath logs/{output[0]}.err)'
#dep_list='"$(sed '\''s/ /,/g'\'' <<< '\''{dependencies}'\'' )"'
cuda=
#cuda='$(test "{resources.cuda}" -eq 1 && echo -l cuda=1)'

#cluster_cmd="qsub -P shendure_fowler -terse -l mfree=${mem_arg}M -l h_rt=48:0:0 -l 'hostname=!s022&!s025' -o $out_path -e $err_path -pe serial {threads}"
#cluster_cmd="echo {log}"
cluster_cmd="qsub -terse -l mfree=${mem_arg}M -l h_rt=48:0:0 -l h=fl004 -o $out_path -e $err_path $cuda -pe serial {threads}"

configfile=
if test -f config.yaml; then
    configfile='--configfile  config.yaml'
fi
if test -t default-config.yaml; then
    configfile='--configfile  default-config.yaml'
fi

snakemake \
    --cluster "$cluster_cmd" \
    --cluster-cancel "qdel" \
    -j 128 \
    $* --cores 50 --resources mem_mb=2000000 --set-resource-scopes mem_mb=global threads=global \
    --rerun-triggers mtime params input software-env \
    --keep-incomplete --latency-wait 120 $configfile \
    --default-resources cuda=0
    #--immediate-submit --notemp \

