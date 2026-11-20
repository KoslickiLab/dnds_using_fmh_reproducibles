#!/usr/bin/env bash
# Copy everything the manuscript uses from this folder into the manuscript folder:
#   results/figures/panels/mtb_*.pdf   ->  <manuscript>/figures/        (one file per panel)
#   results/numbers/mtb.tex            ->  <manuscript>/bioapp/numbers/
#   latex/mtb_*.tex                    ->  <manuscript>/bioapp/         (text, captions, figure layout)
# The six references of latex/mtb.bib are added to the manuscript's bibliography by hand.
#
#   usage:  bash export_to_manuscript.sh <manuscript folder>
set -euo pipefail
manuscript=${1%/}
mkdir -p "$manuscript/figures" "$manuscript/bioapp/numbers"
cp results/figures/panels/mtb_*.pdf "$manuscript/figures/"
cp results/numbers/mtb.tex "$manuscript/bioapp/numbers/"
cp latex/mtb_text.tex latex/mtb_captions.tex latex/mtb_figures.tex "$manuscript/bioapp/"
echo "copied $(ls results/figures/panels/mtb_*.pdf | wc -l) panels, the number file and the vignette's LaTeX to $manuscript"
