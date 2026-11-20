#!/bin/bash
# The published data the investigation starts from, apart from AllTheBacteria (steps 01 and 03).
# Every file is checked against the MD5 sum of the copy used for the vignette.
#
#   usage:  bash analysis/01b_get_published_data.sh [data_folder]      (default: data)
#
#   <data>/gan2025/supp/files/   supplement of Gan et al. 2025 (mSystems 10:e00898-25, PMC12542625), from
#                                Europe PMC: Table S1 (runs and sublineages, s0008), Tables S2 and S3 (PE/PPE
#                                genes and their high-quality segments, s0009 and s0010)
#   <data>/essential/files/      supplement of DeJesus et al. 2017 (mBio 8:e02133-16, PMC5241402): Table S3
#                                (essentiality calls, mbo002173137st3.xlsx)
#   <data>/reference/            the inferred MTBC ancestral sequence (Comas et al. 2013) as distributed with
#                                the code of Gan et al. (github.com/irongump/MTB_PEPPE), and the H37Rv RefSeq
#                                assembly GCF_000195955.2 (annotation and sequence)
set -euo pipefail
DATA=${1:-data}
mkdir -p $DATA/gan2025/supp $DATA/essential $DATA/reference
fetch() { [ -s "$2" ] || curl -sSfL --retry 5 -o "$2" "$1"; }

fetch https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12542625/supplementaryFiles $DATA/gan2025/supp/europepmc_supp.zip
unzip -q -o $DATA/gan2025/supp/europepmc_supp.zip -d $DATA/gan2025/supp/files
fetch https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5241402/supplementaryFiles $DATA/essential/dejesus2017_supp.zip
unzip -q -o $DATA/essential/dejesus2017_supp.zip -d $DATA/essential/files
fetch https://raw.githubusercontent.com/irongump/MTB_PEPPE/4cf810102596cd369f03ac1a6d82d64b2cb5bd01/tb.ancestor.fasta \
      $DATA/reference/tb.ancestor.fasta
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2
for f in GCF_000195955.2_ASM19595v2_genomic.gff.gz GCF_000195955.2_ASM19595v2_genomic.fna.gz; do
    fetch $NCBI/$f $DATA/reference/$f
done

md5sum -c - <<EOF
ee39b42c3767b0e169e65aa1c1fb6465  $DATA/gan2025/supp/files/msystems.00898-25-s0008.xlsx
8e238785ad5dd5f6ec35f974cc67982d  $DATA/gan2025/supp/files/msystems.00898-25-s0009.xlsx
97f5685a893f27d21824027b224e9a44  $DATA/gan2025/supp/files/msystems.00898-25-s0010.xlsx
016df0caf2991dbfe3cb9fae0f90b321  $DATA/essential/files/mbo002173137st3.xlsx
adb84f76313a09a63cddb5b10a0df5ab  $DATA/reference/tb.ancestor.fasta
4d53755c7968f2509d772b12ba2c4136  $DATA/reference/GCF_000195955.2_ASM19595v2_genomic.gff.gz
f93726a0a6e57af2bf1173998fbeb546  $DATA/reference/GCF_000195955.2_ASM19595v2_genomic.fna.gz
EOF
