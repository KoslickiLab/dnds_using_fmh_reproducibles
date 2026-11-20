#!/bin/bash
# Extract the AllTheBacteria 0.2 assemblies of the isolates of Gan et al. 2025 from the 32 downloaded
# batch archives (step 03). Only the samples listed in data/isolates/gan_to_atb.tsv are extracted.
#   out: data/assemblies/<sample>.fa   (one contig set per isolate; about 4.4 MB each)
set -euo pipefail
source "$(dirname "$0")/common.sh"; start_log 10_extract_assemblies
A=$MTB/data/atb/0.2/assembly
OUT=$MTB/data/assemblies; mkdir -p $OUT $MTB/work/members
# sample lists per batch (column 7 is the batch file name, column 4 the ATB sample accession)
awk -F'\t' 'NR>1 && $7!="" && $4!="" {b=$7; sub(/\.asm\.tar\.xz$/,"",b); print b"/"$4".fa" >> ("'"$MTB"'/work/members/" b ".txt")}' \
  $MTB/data/isolates/gan_to_atb.tsv
for m in $MTB/work/members/*.txt; do sort -u $m -o $m; done
echo "$(ls $MTB/work/members/*.txt | wc -l) batches, $(cat $MTB/work/members/*.txt | wc -l) samples requested"
# one decompression per batch, 16 at a time (xz is the bottleneck; 2 threads each)
ls $MTB/work/members/*.txt | xargs -P 16 -I{} bash -c '
  m={}; b=$(basename $m .txt)
  tar -I "xz -T2" -xf '"$A"'/$b.asm.tar.xz -C '"$OUT"' --strip-components=1 -T $m 2>/dev/null || echo "WARN batch $b incomplete"'
echo "assemblies extracted: $(ls $OUT | wc -l)"
du -sh $OUT
