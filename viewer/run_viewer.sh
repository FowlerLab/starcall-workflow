#!/bin/bash
# Starts the dot/cell viewer for a pipeline run.
#
#   interactively, on a compute node from qlogin:
#       viewer/run_viewer.sh [run_dir] [server options]
#   or as a batch job (the tunnel command is printed to logs/viewer.out):
#       qsub -cwd -o logs/viewer.out -j y viewer/run_viewer.sh [run_dir] [server options]
#
# Server options: --well well1 (starting well)  --port 8000  --host 0.0.0.0  --token new (new access token)  --prebuild-cache
# An access token is always required; it is saved in ~/.config/starcall-viewer/token.
# See viewer/README.md (Access).
#$ -l mfree=4G
#$ -l h_rt=12:0:0
#$ -pe serial 2

source /net/fowler/vol1/shared/miniconda3/etc/profile.d/conda.sh
conda activate /net/fowler/vol1/shared/miniconda3/envs/starcall-viewer


# the repo directory (holding viewer/), from this script's location or the qsub working directory
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)"
test -d "$repo_dir/viewer" || repo_dir="$(pwd)"

run_dir="$(pwd)"
if test $# -gt 0 && test "${1#-}" = "$1"; then
    run_dir="$1"
    shift
fi

cd "$repo_dir" && exec env PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1 "/net/fowler/vol1/shared/miniconda3/envs/starcall-viewer/bin/python" -m viewer.server "$run_dir" "$@"
