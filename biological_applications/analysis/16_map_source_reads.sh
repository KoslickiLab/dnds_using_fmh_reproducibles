#!/usr/bin/env bash
# Step 16 (optional, large). Map the metagenome from which Methanobrevibacter sp. RGIG2411 was
# assembled onto the two genomes of the pair, to ask which organism the shared region belongs to.
#
#   usage:  bash analysis/16_map_source_reads.sh pairs/methanobrevibacter_saccharibacteria
#
# The run (settings.ini: source_run) holds 106 million read pairs: the download needs about
# 100 GB of disk in data/reads/ and the mapping about one hour on 32 threads
# (THREADS=32 bash analysis/16_map_source_reads.sh ...). The results of this step are in the
# repository (results/read_depth.tsv, results/read_pairs.tsv), so it only needs to be run to
# reproduce them from the raw reads.
#
# Reference: genome A (the bacterium, in which the shared region lies inside a long scaffold)
# plus genome B (the archaeon) WITHOUT its contigs that consist of shared sequence only
# (results/partner_contigs.tsv, not_shared_bp = 0). The shared region is then present exactly
# once in the reference, so every read from it, whichever organism it comes from, maps there.
#
# Output: work/<pair>/reads/mapped.bam (+ index), then step 17 summarizes it.
set -euo pipefail
pair_dir=${1%/}
threads=${THREADS:-8}
run=$(awk -F' *= *' '$1=="source_run"{print $2}' "$pair_dir/settings.ini")
genome_a=$(awk -F'\t' '$1=="A"{print $2}' "$pair_dir/genomes.tsv")
genome_b=$(awk -F'\t' '$1=="B"{print $2}' "$pair_dir/genomes.tsv")
work=work/$(basename "$pair_dir")/reads
mkdir -p "$work" data/reads

# 1. the reads
if [ ! -s "data/reads/${run}_2.fastq" ]; then
    echo "downloading $run"
    prefetch "$run" --max-size 100G --output-directory data/reads
    fasterq-dump "data/reads/$run" --split-files --threads "$threads" --outdir data/reads --temp "$work"
    rm -r "data/reads/${run:?}"
fi

# 2. the reference: genome A + genome B without its wholly shared contigs
awk -F'\t' 'NR>1 && $4==0 {print $1}' "$pair_dir/results/partner_contigs.tsv" > "$work/wholly_shared_contigs.txt"
echo "contigs of genome B left out of the reference: $(tr '\n' ' ' < "$work/wholly_shared_contigs.txt")"
awk 'NR==FNR {skip[$1]; next}  /^>/ {keep = !(substr($1, 2) in skip)}  keep' \
    "$work/wholly_shared_contigs.txt" "data/genomes/$genome_b/genome.fna" > "$work/genome_b_without_shared_contigs.fna"
cat "data/genomes/$genome_a/genome.fna" "$work/genome_b_without_shared_contigs.fna" > "$work/reference.fna"
bowtie2-build --quiet --threads "$threads" "$work/reference.fna" "$work/reference"

# 3. map the read pairs; reads that do not align are not kept
bowtie2 --threads "$threads" -x "$work/reference" -1 "data/reads/${run}_1.fastq" -2 "data/reads/${run}_2.fastq" \
        --no-unal 2> "$work/bowtie2.log" \
    | samtools sort --threads 8 -m 2G -o "$work/mapped.bam" -
samtools index "$work/mapped.bam"
tail -1 "$work/bowtie2.log"
