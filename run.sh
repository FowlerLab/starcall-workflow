#!/bin/bash

test -e logs/latest && rm -r logs/latest
mkdir logs/latest

mem_arg='$(expr {resources.mem_mb} / {threads})'
#out_path=$(realpath ./logs/latest/)
#err_path=$(realpath ./logs/latest/)
#out_path='$(realpath {log})'
#err_path='$(realpath {log})'
out_path='$(mkdir -p $(dirname {log}.out); realpath {log}.out)'
err_path='$(mkdir -p $(dirname {log}.err); realpath {log}.err)'
#dep_list='"$(sed '\''s/ /,/g'\'' <<< '\''{dependencies}'\'' )"'

#cluster_cmd="qsub -P shendure_fowler -terse -l mfree=${mem_arg}M -l h_rt=48:0:0 -l 'hostname=!s022&!s025' -o $out_path -e $err_path -pe serial {threads}"
#cluster_cmd="echo {log}"
cluster_cmd="qsub -terse -l mfree=${mem_arg}M -l h_rt=48:0:0 -o $out_path -e $err_path -pe serial {threads}"

snakemake \
    --cluster "$cluster_cmd" \
    --cluster-cancel "qdel" \
    -j 128 \
    $* --cores 50 --resources mem_mb=200000 --set-resource-scopes mem_mb=global threads=global \
    --rerun-triggers mtime params input software-env \
    --keep-incomplete --latency-wait 30
    #--immediate-submit --notemp \

