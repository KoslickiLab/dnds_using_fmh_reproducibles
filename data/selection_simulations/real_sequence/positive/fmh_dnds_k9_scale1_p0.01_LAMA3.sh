#!/bin/bash
set -eoux pipefail

#Parameters
wd=/data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/positive/fmh_k9
k='9'
scaled=1
translate_cds='no' #indicate with yes or no
m="sngl" #branchwater mode, single fasta file use sourmash sketch, multiple fasta files to use sourmash sketch
c=1

#Note that if using branchwater mode, makesure working directory is at the beginning of filename. THIS IS NOT REQUIERED FOR OTHER MODES
#cds_input_list=${wd}/dataset.csv

#using sngl mode
dna_filename=${wd}/../positive_selection_queries_0.01_K06240.fna
protein_filename=${wd}/../positive_selection_translated_queries_0.01_K06240.faa

#outname
out=K06240_positive_0.01

#Run FMH Omega
python3 /data/jzr5814/repositories/dnds-using-fmh/src/script_fmh_omega.py --dna_fasta ${dna_filename} --protein_fasta ${protein_filename} --scaled_input ${scaled} --ksize ${k} --mode ${m} --translate ${translate_cds} --outname ${out} --directory ${wd} --cores ${c} 



