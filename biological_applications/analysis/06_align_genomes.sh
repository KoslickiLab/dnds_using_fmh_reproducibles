#!/usr/bin/env bash
# Step 6. Whole-genome alignment (nucmer, default settings) of genome A of the pair against
# genome B and against every other genome listed in the pair's genomes.tsv.
#
#   usage:  bash analysis/06_align_genomes.sh <pair folder>
#
# Output: work/<pair folder name>/alignments/<accession>.coords
#         one line per aligned block: positions in A, positions in the other genome, lengths,
#         percent identity, contig lengths, coverages, contig names (show-coords -THrcl)
set -euo pipefail
pair_dir=${1%/}
out=work/$(basename "$pair_dir")/alignments
mkdir -p "$out"

genome_a=$(awk -F'\t' '$1=="A"{print $2}' "$pair_dir/genomes.tsv")

for accession in $(awk -F'\t' 'NR>1 && $1!="A"{print $2}' "$pair_dir/genomes.tsv" | sort -u); do
    if [ -s "$out/$accession.coords" ]; then continue; fi
    echo "aligning $genome_a (reference) with $accession"
    nucmer --threads "${THREADS:-8}" --prefix "$out/$accession" \
        "data/genomes/$genome_a/genome.fna" "data/genomes/$accession/genome.fna" 2> "$out/$accession.nucmer.log"
    show-coords -THrcl "$out/$accession.delta" > "$out/$accession.coords"
done
