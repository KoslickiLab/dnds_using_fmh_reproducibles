#!/bin/bash
# Runs analysis/07_reads_to_snvs.sh for every run in data/isolates/read_pipeline_runs.txt
# (300 calibration runs that also have an ATB assembly, then the 10,570 runs without one; step 06).
# 80 runs at a time with 3 threads each (240 threads). Rerunnable: finished runs are skipped.
set -euo pipefail
source "$(dirname "$0")/common.sh"; start_log 08_run_read_pipeline
cat $MTB/data/isolates/read_pipeline_runs.txt | xargs -P ${JOBS:-80} -I{} sh -c "$MTB/analysis/07_reads_to_snvs.sh {} 3 || echo 'FAILED {}'"
echo "done: $(ls $MTB/work/snv_reads/*.stats.tsv | wc -l) runs with stats"
