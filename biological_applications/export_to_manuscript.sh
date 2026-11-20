#!/usr/bin/env bash
# Copy everything the manuscript uses from this repository into the manuscript folder:
#   output/figures/*.pdf   ->  <manuscript>/figures/
#   output/numbers/*.tex   ->  <manuscript>/bioapp/numbers/
#   output/data/*          ->  <manuscript>/supplementary_data/
#
#   usage:  bash export_to_manuscript.sh <manuscript folder>
set -euo pipefail
manuscript=${1%/}
mkdir -p "$manuscript/figures" "$manuscript/bioapp/numbers" "$manuscript/supplementary_data"
cp output/figures/*.pdf "$manuscript/figures/"
cp output/numbers/*.tex "$manuscript/bioapp/numbers/"
cp output/data/* "$manuscript/supplementary_data/"
echo "copied $(ls output/figures/*.pdf | wc -l) figures, $(ls output/numbers/*.tex | wc -l) number files and $(ls output/data | wc -l) data files to $manuscript"
