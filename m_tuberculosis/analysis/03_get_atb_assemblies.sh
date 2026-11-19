#!/bin/bash
# AllTheBacteria 0.2 assembly batches (tar.xz) that hold the isolates of Gan et al. 2025, checked against
# the release's MD5 sums. Step 10 extracts the isolates' assemblies from them.
set -euo pipefail
source "$(dirname "$0")/common.sh"; start_log 03_get_atb_assemblies
B=https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/assembly
D=$MTB/data/atb/0.2/assembly; mkdir -p $D $MTB/data/assemblies
cd $D
files=$(cut -f7 $MTB/data/isolates/gan_to_atb.tsv | tail -n +2 | grep -v '^$' | sort -u)
echo "$(echo "$files" | wc -l) batches needed"
echo "$files" | xargs -P 8 -I{} sh -c '[ -s {} ] || curl -sSfL -o {} '$B'/{}'
grep -F "$(echo "$files" | sed 's#^#./assembly/#')" ../md5sum.txt | sed 's#  \./assembly/#  #' | md5sum -c - | grep -vc ': OK$' | xargs -I{} echo "{} batches failed md5"
# the assemblies of the isolates are extracted from these archives in step 10
du -sh $MTB/data/atb/0.2/assembly
