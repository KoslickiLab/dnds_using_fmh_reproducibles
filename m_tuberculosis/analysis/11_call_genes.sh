#!/bin/bash
# Prodigal gene calls for every extracted assembly, with the settings used for the GTDB genomes of the
# main manuscript (prodigal -m -p single, translation table 11). The FracMinHash dN/dS estimator needs
# the coding sequences and their translations, nothing else: no functional annotation is done anywhere
# in this investigation.
#   out: work/genes/<sample>.ffn.gz   coding sequences
#        work/genes/<sample>.faa.gz   their translations
#        work/genes/<sample>.log      prodigal stderr (kept only if prodigal failed)
# Rerunnable: samples that already have both outputs are skipped.
set -euo pipefail
source "$(dirname "$0")/common.sh"; start_log 11_call_genes
OUT=$MTB/work/genes; mkdir -p $OUT
export OUT
one() {
  s=$(basename "$1" .fa)
  [ -s "$OUT/$s.faa.gz" ] && [ -s "$OUT/$s.ffn.gz" ] && return 0
  prodigal -i "$1" -m -p single -g 11 -q -f gff -o /dev/null \
      -d "$OUT/$s.ffn.tmp" -a "$OUT/$s.faa.tmp" 2> "$OUT/$s.log" || { echo "FAILED $s"; return 0; }
  gzip -1 -c "$OUT/$s.ffn.tmp" > "$OUT/$s.ffn.gz" && gzip -1 -c "$OUT/$s.faa.tmp" > "$OUT/$s.faa.gz"
  rm -f "$OUT/$s.ffn.tmp" "$OUT/$s.faa.tmp" "$OUT/$s.log"
}
export -f one
find $MTB/data/assemblies -name '*.fa' | xargs -P ${JOBS:-280} -I{} bash -c 'one {}'
echo "gene calls: $(find $OUT -name '*.faa.gz' | wc -l) of $(find $MTB/data/assemblies -name '*.fa' | wc -l) assemblies"
du -sh $OUT
