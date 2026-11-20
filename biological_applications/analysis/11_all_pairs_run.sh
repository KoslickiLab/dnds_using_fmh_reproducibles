#!/usr/bin/env bash
# Step 11. Steps 1-5 for every genome pair of the hierarchical edge bundling figure
# (all_pairs/pairs.tsv: 140 pairs, 243 genomes).
#
#   usage:  bash analysis/11_all_pairs_run.sh
#
#   JOBS=4 THREADS=2 bash analysis/11_all_pairs_run.sh     # pairs run in parallel, threads per pair
#
# Each pair gets a folder work/all_pairs/<pair>/ with a two-line genomes.tsv, and the same
# scripts that analyse the pairs examined in detail are run on it:
#   step 3  near-identical genes            -> results/genes_A.tsv, genes_B.tsv
#   step 4  FracMinHash dN/dS               -> results/fmh_dnds.tsv
#   step 5  alignment-based dN/dS           -> results/orthologs.tsv
# Step 12 then condenses these into all_pairs/results/all_pairs_summary.tsv.
# A pair that is already finished is skipped, so the script can be re-run after an interruption.
# Expect about 0.7 GB of downloads and one to two CPU-minutes per pair (each pair needs up to 2 GB of memory).
set -euo pipefail
export JOBS=${JOBS:-4}
export THREADS=${THREADS:-2}
export PYTHONPATH=analysis

# genomes are shared between pairs, so they are downloaded and annotated once, before the pairs are run
mkdir -p work/all_pairs
{ printf 'role\taccession\n'; tail -n +2 all_pairs/genomes.tsv | cut -f1 | sed 's/^/genome\t/'; } > work/all_pairs/genomes.tsv
bash analysis/01_download_genomes.sh work/all_pairs
THREADS=$((JOBS * THREADS)) python analysis/02_call_genes.py work/all_pairs

run_one_pair() {
    pair=$1; genome_1=$2; genome_2=$3
    dir=work/all_pairs/$pair
    if [ -s "$dir/results/orthologs.tsv" ]; then return; fi
    mkdir -p "$dir"
    printf 'role\taccession\nA\t%s\nB\t%s\n' "$genome_1" "$genome_2" > "$dir/genomes.tsv"
    {
        python analysis/03_shared_genes.py "$dir"
        python analysis/04_fmh_dnds.py "$dir"
        python analysis/05_conventional_dnds.py "$dir"
    } > "$dir/log.txt" 2>&1 || echo "FAILED: $pair (see $dir/log.txt)"
    echo "finished $pair"
}
export -f run_one_pair
tail -n +2 all_pairs/pairs.tsv | cut -f1-3 | tr '\t' ' ' | xargs -P "$JOBS" -L 1 bash -c 'run_one_pair "$@"' _
echo "pairs finished: $(ls work/all_pairs/*/results/orthologs.tsv | wc -l) of $(($(wc -l < all_pairs/pairs.tsv) - 1))"
