#!/usr/bin/env bash
# Reproduce every number, table and figure of the M. tuberculosis vignette from public data.
#
#   conda env create -f environment.yml && conda activate fmh_dnds_mtb
#   bash reproduce_all.sh
#
# This is a large computation, meant for a server: it downloads 16 GB of AllTheBacteria assembly
# archives (165 GB once the 39,379 assemblies are extracted) and needs about 5 TB of disk for work/
# (the sketches at scale factor 1 and the per-pair tables of 775 million pairs per group of genes).
# Step 32 ran on 250 cores; the FracMinHash dN/dS tool's own run times are in
# results/trends_tool/timings.tsv and results/tool/sketching.tsv. Everything is written inside this
# directory. The optional read pipeline (steps 06-08 and 34) is not run; see README.md.
set -euo pipefail
cd "$(dirname "$0")"
export JOBS=${JOBS:-200}      # processes for the parallel steps
export CORES=${CORES:-250}    # cores for the FracMinHash dN/dS runs (step 32)
export THREADS=${THREADS:-48} # cores for the draft-assembly runs of step 12
export PYTHONDONTWRITEBYTECODE=1

# ---- data
bash   analysis/00_get_tools.sh
bash   analysis/01_get_atb_metadata.sh
bash   analysis/01b_get_published_data.sh
python analysis/02_match_gan_to_atb.py
bash   analysis/03_get_atb_assemblies.sh
python analysis/05_ena_runinfo.py
bash   analysis/10_extract_assemblies.sh

# ---- sequences of every isolate
bash   analysis/11_call_genes.sh                 # gene calls of the draft assemblies (for step 12 only)
python analysis/12_smoke_scales.py 2             # draft assemblies compared directly, main mode
python analysis/13_gene_groups.py
bash   analysis/14_variants_from_assemblies.sh
python analysis/14b_coverage.py
python analysis/24_site_ratio.py
python analysis/30_write_inputs.py --jobs $JOBS

# ---- FracMinHash dN/dS of every pair, after checking reference mode against main mode
python checks/31a_check_inputs.py
bash   checks/31_check_reference_mode.sh
python checks/31b_check_summary.py
bash   analysis/32_run_reference_mode.sh
python analysis/33_conventional.py
python analysis/35_sketch_times.py

# ---- summaries, figures, numbers
python analysis/40_summarise_tool_output.py --jobs $JOBS
python analysis/40b_exact_quantiles.py --jobs $JOBS
python analysis/41_trends.py
python analysis/44_divergence_tail.py
python figures/42_figures.py
python analysis/43_numbers.py

# ---- every panel recomputed from the raw output, independently of the steps above
python checks/46_review_figures.py --jobs $JOBS
