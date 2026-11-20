#!/bin/bash
# FracMinHash dN/dS for every pair of isolates, with the lab's tool (external/dnds-using-fmh,
# fmh_dnds.py --mode reference), from the FASTA files of step 30.
#
#   1. whole coding genome, every isolate (assembly- and read-derived): the tool sketches all of them
#      at scale factor 1 (every k-mer), in 250 shards side by side (--sketch_shards), and compares
#      every pair.
#   2. whole coding genome, the assembly-derived isolates only, from the same sketches, at scale
#      factor 1 and downsampled to 2, 3, 5, 10, 100 and 1000.
#   3. each of the four gene groups, assembly-derived isolates: sketched at scale factor 1 and
#      compared at 1, 2, 3, 5, 10, 100 and 1000.
# Parameters are the manuscript's: k_aa = 7 (k_nt = 21), max containment, constant 0.77/2.23.
# --threshold 0 so that every pair is reported (all containments here are above 0.9 anyway).
#
# Rerunnable: a run whose shard sketches are all there uses them (--sketches) instead of sketching
# again. How long sketching took is read from this log by step 35.
#
#   out: work/tool/<run>/sketch_shards/<shard>/data.zip     the sketches
#        work/tool/<run>/reference_scaled<S>/pairs/*.parquet   one row per pair, main mode's columns
#        work/tool/<run>/reference_scaled<S>/run.json, names.txt
set -euo pipefail
source "$(dirname "$0")/common.sh"; start_log 32_run_reference_mode
TOOL=$MTB/external/dnds-using-fmh/fmh_dnds.py
export PYTHONDONTWRITEBYTECODE=1
CORES=${CORES:-250}
SHARDS=${SHARDS:-250}
SWEEP=2,3,5,10,100,1000
common="--mode reference --reference ancestral_reference --ksize 7 --threshold 0 --cores $CORES"
cd $(dirname $TOOL)
mkdir -p $MTB/work/tool/{all_cds_every_isolate,all_cds,essential,nonessential,pe_ppe,pe_ppe_hq}

# the sketches of a run: the existing shard zips if all of them are there, otherwise sketch the csv
sketches_of() {
  local run=$1 genomes=$2
  local zips=($MTB/work/tool/$run/sketch_shards/*/data.zip)
  if [ -e "${zips[0]}" ] && [ ${#zips[@]} -eq $SHARDS ]; then
    echo "--sketches ${zips[*]} --scaled_input 1"
  else
    echo "--sketch_shards $SHARDS --scaled_input 1 --fasta_input_list $genomes"
  fi
}

echo "=== 1. whole coding genome, every isolate $(date -Is)"
python $TOOL $common $(sketches_of all_cds_every_isolate $MTB/work/fasta/all_cds/genomes.csv) \
    --directory $MTB/work/tool/all_cds_every_isolate

echo "=== 2. whole coding genome, assembly-derived isolates, scale factors 1 and $SWEEP $(date -Is)"
python $TOOL $common --sketches $MTB/work/tool/all_cds_every_isolate/sketch_shards/*/data.zip \
    --names $MTB/work/fasta/assembly_names.txt --scaled_input 1 --downsample_to $SWEEP \
    --directory $MTB/work/tool/all_cds

for group in essential nonessential pe_ppe pe_ppe_hq; do
  echo "=== 3. $group, scale factors 1 and $SWEEP $(date -Is)"
  python $TOOL $common $(sketches_of $group $MTB/work/fasta/$group/genomes.csv) \
      --downsample_to $SWEEP --directory $MTB/work/tool/$group
done
echo "=== done $(date -Is)"
du -sh $MTB/work/tool/*
