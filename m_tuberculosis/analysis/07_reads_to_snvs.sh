#!/bin/bash
# One run: SRA (AWS open data) -> reads (downsampled to at most 100x) -> bwa mem on the MTBC ancestral reference -> fixed SNVs and indels.
# Thresholds follow Gan et al. 2025 (VarScan: mapping quality > 30, depth > 10, alt frequency > 0.9; base quality 20).
#   usage: analysis/07_reads_to_snvs.sh RUN [threads]
#   out:   work/snv_reads/RUN.calls.tsv.gz   pos ref alt depth alt_fraction   (fixed variants, ancestral-reference coordinates)
#          work/snv_reads/RUN.stats.tsv      raw bases, downsampling fraction, mean depth and fraction of genome with depth >= 10 (MAPQ > 30)
set -euo pipefail
source "$(dirname "$0")/env.sh"
RUN=$1; T=${2:-3}
OUT=$MTB/work/snv_reads; mkdir -p $OUT
[ -s $OUT/$RUN.stats.tsv ] && exit 0
W=$MTB/work/read_tmp/$RUN; mkdir -p $W; trap 'rm -rf $W' EXIT
REF=$MTB/data/reference/tb.ancestor.fasta
curl -sSfL --retry 5 --retry-delay 20 -o $W/$RUN.sra "https://sra-pub-run-odp.s3.amazonaws.com/sra/$RUN/$RUN"
fasterq-dump --split-3 --skip-technical -e $T -t $W/fq_tmp -O $W $W/$RUN.sra > /dev/null 2>&1
rm -f $W/$RUN.sra
reads=$(ls $W/${RUN}_1.fastq $W/${RUN}_2.fastq 2>/dev/null || ls $W/$RUN.fastq)
# downsample runs deeper than 100x (same seed for both mates keeps the pairs)
bases=$(cat $reads | awk 'NR%4==2{s+=length($0)} END{print s}')
frac=$(awk -v b=$bases 'BEGIN{f=100*4411532/b; print (f<1 ? f : 1)}')
if awk -v f=$frac 'BEGIN{exit !(f<1)}'; then
  for r in $reads; do seqtk sample -s 11 $r $frac > $r.sub && mv $r.sub $r; done
fi
bwa mem -t $T $REF $reads 2>/dev/null \
  | samtools sort -@ 1 -m 1500M -T $W/sort -o $W/$RUN.bam - 2>/dev/null
rm -f $W/*.fastq
samtools index $W/$RUN.bam
# pileup of the whole genome; depth capped at 500 reads per position
bcftools mpileup -f $REF -q 31 -Q 20 -d 500 -a FORMAT/AD,FORMAT/DP --threads 1 $W/$RUN.bam 2>/dev/null \
  | bcftools call -m --ploidy 1 -v -Ou 2>/dev/null \
  | bcftools query -f '%POS\t%REF\t%ALT\t[%DP]\t[%AD]\n' \
  | awk -F'\t' 'BEGIN{OFS="\t"} {n=split($5,a,","); s=0; for(i=1;i<=n;i++) s+=a[i]; if (s==0) next; split($3,alt,",");
        best=2; for(i=3;i<=n;i++) if (a[i]>a[best]) best=i; f=a[best]/s;
        if ($4>10 && f>0.9) print $1,$2,alt[best-1],$4,sprintf("%.3f",f)}' | gzip > $OUT/$RUN.calls.tsv.gz
samtools depth -a -Q 31 $W/$RUN.bam | awk -v run=$RUN -v b=$bases -v f=$frac 'BEGIN{OFS="\t"} {n++; s+=$3; if ($3>=10) c++} END{print "run","raw_bases","downsampled_to","mean_depth","fraction_ge10"; print run, b, f, s/n, c/n}' > $OUT/$RUN.stats.tmp
mv $OUT/$RUN.stats.tmp $OUT/$RUN.stats.tsv
