#!/bin/bash
# Substitutions of every assembly against the inferred MTBC ancestral sequence (Comas et al. 2013),
# which has the coordinates of H37Rv (NC_000962.3, 4,411,532 bp) and is the reference Gan et al. 2025
# called their variants against.
#
# Why this step exists: comparing the gene calls of two independent draft assemblies measures their
# difference in gene content as much as their difference in sequence (see results/smoke/summary.tsv),
# and that difference removes nucleotide and protein k-mers in equal proportion, which is exactly the
# regime in which the estimator's denominator collapses. Anchoring every isolate to one gene set
# removes that term. Only substitutions are kept, so every isolate has the same coding coordinates.
#
#   out: work/variants/<sample>.tsv.gz   pos(1-based) ref alt, for substitutions with mapping quality >= 20
#        work/variants/<sample>.cov      covered bases of the reference
# Rerunnable: samples that already have both outputs are skipped.
set -euo pipefail
source "$(dirname "$0")/common.sh"; start_log 14_variants_from_assemblies
REF=$MTB/data/reference/tb.ancestor.fasta
OUT=$MTB/work/variants; mkdir -p $OUT
export OUT REF
one() {
  s=$(basename "$1" .fa)
  [ -s "$OUT/$s.cov" ] && [ -s "$OUT/$s.tsv.gz" ] && return 0
  minimap2 -cx asm5 --cs -t 2 "$REF" "$1" 2>/dev/null | sort -k6,6 -k8,8n \
    | paftools.js call -l 2000 -L 5000 - 2>/dev/null > "$OUT/$s.call" || { echo "FAILED $s"; return 0; }
  awk -F'\t' 'BEGIN{OFS="\t"}
      $1=="V" && $6>=20 && length($7)==1 && length($8)==1 && $7~/^[acgtACGT]$/ && $8~/^[acgtACGT]$/ {
        print $3+1, toupper($7), toupper($8) }' "$OUT/$s.call" | gzip -1 -c > "$OUT/$s.tsv.gz"
  awk -F'\t' '$1=="R"{c+=$4-$3} END{printf "%d\n", c+0}' "$OUT/$s.call" > "$OUT/$s.cov"
  rm -f "$OUT/$s.call"
}
export -f one
find $MTB/data/assemblies -name '*.fa' | xargs -P ${JOBS:-140} -I{} bash -c 'one {}'
echo "isolates with variants: $(find $OUT -name '*.cov' | wc -l)"
find $OUT -name '*.cov' -exec cat {} + | awk '{s+=$1; n++} END{if(n) print "mean covered bases:", s/n, "of 4411532 =", 100*s/n/4411532 "%"}'
du -sh $OUT
