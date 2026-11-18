#!/usr/bin/env bash
# Run every analysis step for one genome pair examined in detail, in order.
#
#   usage:  bash analysis/run_pair.sh pairs/<name>
#
# Run from the folder that contains README.md, with the conda environment active
# (conda activate fmh_dnds_applications) and after `bash analysis/00_get_lab_code.sh`.
# Each step reads the files written by the steps before it; see the header of each script.
set -euo pipefail
pair_dir=${1%/}
export PYTHONPATH=analysis

bash   analysis/01_download_genomes.sh            "$pair_dir"   # genomes and annotation from NCBI
python analysis/02_call_genes.py                  "$pair_dir"   # Prodigal gene calls (as in GTDB)
python analysis/03_shared_genes.py                "$pair_dir"   # which genes are near-identical between the two genomes
python analysis/04_fmh_dnds.py                    "$pair_dir"   # FracMinHash dN/dS with and without those genes
python analysis/05_conventional_dnds.py           "$pair_dir"   # alignment-based dN/dS of the orthologs
bash   analysis/06_align_genomes.sh               "$pair_dir"   # nucmer: genome A against genome B and the other genomes
python analysis/07_shared_regions.py              "$pair_dir"   # shared regions and their annotated genes
python analysis/08_scaffold_gaps.py               "$pair_dir"   # scaffold gaps next to the shared regions
python analysis/09_presence_in_other_genomes.py   "$pair_dir"   # presence of each region in the other genomes
python analysis/10_pair_numbers.py                "$pair_dir"   # output/numbers/<name>.tex
