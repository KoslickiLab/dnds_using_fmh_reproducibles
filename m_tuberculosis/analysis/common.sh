# Settings shared by the shell steps. Source, do not execute.
# Every step runs from the folder that holds README.md, with the conda environment of environment.yml
# active (environment_reads.yml for the read pipeline, steps 07 and 08).
export MTB=$(pwd)
export TMPDIR=${TMPDIR:-$MTB/tmp}
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$TMPDIR" "$MTB/logs"
# write everything a step prints to logs/<step>.log as well
start_log() { exec > >(tee -a "$MTB/logs/$1.log") 2>&1; echo "=== $1 start $(date -Is) on $(hostname) ==="; }
