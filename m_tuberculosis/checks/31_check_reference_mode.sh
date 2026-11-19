#!/bin/bash
# Faithfulness of the tool's reference mode on this collection, with the lab's own check
# (external/dnds-using-fmh/tests/check_reference_mode.py): for
# isolates chosen at random (seed 20261118) plus the reference, main mode (unchanged) and reference
# mode must report identical values for every pair and every column.
#   pe_ppe_hq: 1,000 isolates. all_cds: the first 80 of the same 1,000, because main mode cannot
#   read back its own sketch file once it exceeds 4 GB (manysketch then writes a ZIP64 archive that
#   sourmash scripts pairwise does not recognise); 81 whole coding genomes at scale factor 1 are 3.5 GB.
#   * Reference mode makes its own sketches, in 20 shards side by side (--sketch_shards), so the
#     sharded sketching used in production is checked against main mode's single sketching run.
#   * At each extra scale factor, main mode sketches the genomes again at that scale factor, while
#     reference mode downsamples its scale-factor-1 sketches.
#   pe_ppe_hq (fewest k-mers, fewest mutations): 1,000 isolates, scale factors 1, 2 and 5
#   all_cds   (most k-mers):                     80 isolates, scale factors 1 and 5
# Main-mode runs already finished are reused (--reuse_main); main mode does not depend on the change.
#   usage: bash checks/31_check_reference_mode.sh [group ...]        default: pe_ppe_hq all_cds
#   out:   work/equivalence/<group>/   logs/31_check_reference_mode.log
set -euo pipefail
source "$(dirname "$0")/../analysis/common.sh"; start_log 31_check_reference_mode
TOOL=$MTB/external/dnds-using-fmh
export PYTHONDONTWRITEBYTECODE=1
for group in ${@:-pe_ppe_hq all_cds}; do
  case $group in
    pe_ppe_hq) scales=2,5; isolates=1000 ;;
    all_cds)   scales=5;   isolates=80 ;;
  esac
  echo "--- $group: $isolates isolates and the reference; scale factors 1 and $scales; code $(git -C $TOOL rev-parse --short HEAD)"
  rm -rf $MTB/work/equivalence/$group/reference
  python $TOOL/tests/check_reference_mode.py --fasta_input_list $MTB/work/equivalence/${group}_${isolates}.csv \
      --reference ancestral_reference --directory $MTB/work/equivalence/$group \
      --scaled 1 --downsample_to $scales --cores ${CORES:-200} --threshold 0 \
      --sketch_shards 20 --reuse_main
done
