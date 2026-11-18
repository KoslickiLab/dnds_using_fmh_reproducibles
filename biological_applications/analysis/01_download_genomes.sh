#!/usr/bin/env bash
# Step 1. Download every genome listed in a pair's genomes.tsv from NCBI
# (assembly accessions are version-pinned, so the sequences never change).
#
#   usage:  bash analysis/01_download_genomes.sh <pair folder>
#
# For each accession this creates
#   data/genomes/<accession>/genome.fna       the assembled sequence
#   data/genomes/<accession>/annotation.gff   NCBI's gene annotation (used only for product names)
# Genomes that are already there are skipped, so the step can be re-run safely.
# If the environment variable NCBI_API_KEY is set, the NCBI datasets tool uses it.
set -euo pipefail
pair_dir=$1

# column 2 of genomes.tsv is the accession; the first line is the header
for accession in $(tail -n +2 "$pair_dir/genomes.tsv" | cut -f2 | sort -u); do
    out=data/genomes/$accession
    if [ -s "$out/genome.fna" ]; then continue; fi
    echo "downloading $accession"
    mkdir -p "$out"
    datasets download genome accession "$accession" --include genome,gff3 \
        --filename "$out/ncbi.zip" --no-progressbar
    unzip -q -o "$out/ncbi.zip" -d "$out/unzipped"
    cat "$out"/unzipped/ncbi_dataset/data/"$accession"/*_genomic.fna > "$out/genome.fna"
    # assemblies without an NCBI annotation (many metagenome-assembled genomes) have no gff file
    cp "$out"/unzipped/ncbi_dataset/data/"$accession"/genomic.gff "$out/annotation.gff" 2>/dev/null || echo "  (no NCBI annotation for $accession)"
    rm -r "${out:?}/unzipped" "${out:?}/ncbi.zip"
    sleep 1   # be polite to NCBI
done
