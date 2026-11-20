#!/bin/bash
# AllTheBacteria release 0.2 metadata (assemblies of all ENA/SRA Illumina isolates as of 2023-06-16).
set -euo pipefail
source "$(dirname "$0")/common.sh"; start_log 01_get_atb_metadata
B=https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2
D=$MTB/data/atb/0.2; mkdir -p $D/metadata
cd $D
curl -sSfL -o md5sum.txt $B/md5sum.txt
for f in sample2species2file.tsv.gz ena_metadata.tsv.gz assembly-stats.tsv.gz checkm2.tsv.gz hq_set.sample_list.txt.gz species_calls.tsv.gz sylph.tsv.gz; do
  [ -s metadata/$f ] || curl -sSfL -o metadata/$f $B/metadata/$f
done
grep -E "  (\./)?metadata/" md5sum.txt | grep -E "$(ls metadata | paste -sd'|')" | sed 's#  \./#  #' | md5sum -c -
