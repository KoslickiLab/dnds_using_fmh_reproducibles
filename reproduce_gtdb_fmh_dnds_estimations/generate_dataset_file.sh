#!/bin/bash
set -euo pipefail

# Temporary files
tmp_fna=$(mktemp)
tmp_faa=$(mktemp)

# Find all genome (.fna) files with absolute paths
find . -type f -name "*_protein.fna.gz" -exec realpath {} \; | while read -r file; do
    sample=$(basename "$file" "_protein.fna.gz")
    echo "${sample},${file}"
done | sort > "$tmp_fna"

# Find all protein (.faa) files with absolute paths
find . -type f -name "*_protein.faa.gz" -exec realpath {} \; | while read -r file; do
    sample=$(basename "$file" "_protein.faa.gz")
    echo "${sample},${file}"
done | sort > "$tmp_faa"

# Create the CSV
{
    echo "name,genome_filename,protein_filename"

    join -t, -1 1 -2 1 "$tmp_fna" "$tmp_faa"
} > dataset.csv

# Clean up
