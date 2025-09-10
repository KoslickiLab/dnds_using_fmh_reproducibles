#!/bin/bash
set -eoux pipefail

#Parameters
wd=/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/negative
k='7'
scaled=1
translate_cds='no' #indicate with yes or no
m="sngl" #branchwater mode, single fasta file use sourmash sketch, multiple fasta files to use sourmash sketch
c=1000
t=0

#Note that if using branchwater mode, makesure working directory is at the beginning of filename. THIS IS NOT REQUIERED FOR OTHER MODES
dna=/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/negative_selection_queries_20001_0.01.fna
protein=/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/negative_selection_translated_queries_20001_0.01.faa
out="prate_0.01"

#Run FMH Omega
python3 /data/jzr5814/repositories/dnds-using-fmh/src/script_fmh_omega.py --dna_fasta ${dna} --protein_fasta ${protein} --scaled_input ${scaled} --ksize ${k} --mode ${m} --translate ${translate_cds} --outname ${out} --directory ${wd} --cores ${c} --threshold ${t} 
