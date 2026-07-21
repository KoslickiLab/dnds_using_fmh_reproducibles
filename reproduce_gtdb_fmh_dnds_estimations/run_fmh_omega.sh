#!/bin/bash
set -euo pipefail

source "$1"

python3 /data/jzr5814/repositories/dnds-using-fmh/src/script_fmh_omega.py \
    --fasta_input_list "${wd}/dataset.csv" \
    --scaled_input "$scaled" \
    --ksize "$k" \
    --mode "$mode" \
    --translate "$translate_cds" \
    --outname "$out" \
    --directory "$wd" \
    --cores "$cores" \
    --threshold "$threshold"

