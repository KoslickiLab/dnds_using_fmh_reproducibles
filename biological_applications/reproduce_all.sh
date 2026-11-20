#!/usr/bin/env bash
# Reproduce every number, table and figure of the biological applications.
#
#   conda env create -f environment.yml && conda activate fmh_dnds_applications
#   bash reproduce_all.sh
#
# Needs internet access (about 280 genomes from NCBI, 1 GB) and up to 2 GB of memory per pair
# that runs at the same time. With JOBS=40 the whole run took 20 minutes on our server; with the
# defaults below (4 pairs at a time, 2 threads each) it should take about an hour.
# Everything is written inside this directory. The optional read-level steps 16 and 17 are not
# run (see README.md).
set -euo pipefail
cd "$(dirname "$0")"
export JOBS=${JOBS:-4}          # genome pairs analysed at the same time (step 11)
export THREADS=${THREADS:-2}    # threads per pair
export PYTHONPATH=analysis:figures

bash analysis/00_get_lab_code.sh

# ---- the 140 genome pairs of the hierarchical edge bundling figure (steps 1-5 for each pair)
bash   analysis/11_all_pairs_run.sh
python analysis/12_all_pairs_summary.py

# ---- the pairs examined in detail (steps 1-10 for each pair)
for pair in $(tail -n +2 pairs/examined_pairs.tsv | cut -f2); do
    echo "=================== $pair"
    bash analysis/run_pair.sh "$pair"
done

# ---- experiments and summary numbers
python analysis/13_gene_removal_experiment.py
python analysis/14_translation_table_check.py
python analysis/15_overview_numbers.py
python analysis/18_supplementary_data.py

# ---- figures
bash figures/make_all_figures.sh

# ---- checks: same gene calls and FracMinHash dN/dS values as the manuscript's data
python checks/check_against_gtdb.py
