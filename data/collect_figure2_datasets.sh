#!/bin/bash
set -euo pipefail

mkdir -p figure2_datasets

cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/negative/k5/fmh_omega_5.csv \
   figure2_datasets/fmh_omega_5_negative_5001bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/ksizes/5/negative/fmh_omega_5.csv \
   figure2_datasets/fmh_omega_5_negative_10002bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/negative/k5/fmh_omega_5.csv \
   figure2_datasets/fmh_omega_5_negative_20001bp.csv

cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/negative/fmh_omega_7.csv \
   figure2_datasets/fmh_omega_7_negative_5001bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/fmh_dnds_sketch_protein/negative_selection_redo_sketch_protein_using_faa/dnds.k7.approximations_included.csv \
   figure2_datasets/fmh_omega_7_negative_10002bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/negative/fmh_omega_7.csv \
   figure2_datasets/fmh_omega_7_negative_20001bp.csv

cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/negative/k15/fmh_omega_15.csv \
   figure2_datasets/fmh_omega_15_negative_5001bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/ksizes/15/negative/fmh_omega_15.csv \
   figure2_datasets/fmh_omega_15_negative_10002bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/negative/k15/fmh_omega_15.csv \
   figure2_datasets/fmh_omega_15_negative_20001bp.csv

cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/negative/k21/fmh_omega_21.csv \
   figure2_datasets/fmh_omega_21_negative_5001bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/ksizes/21/negative/fmh_omega_21.csv \
   figure2_datasets/fmh_omega_21_negative_10002bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/negative/k21/fmh_omega_21.csv \
   figure2_datasets/fmh_omega_21_negative_20001bp.csv

cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/positive/k5/fmh_omega_5.csv \
   figure2_datasets/fmh_omega_5_positive_5001bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/ksizes/5/positive/fmh_omega_5.csv \
   figure2_datasets/fmh_omega_5_positive_10002bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/positive/k5/fmh_omega_5.csv \
   figure2_datasets/fmh_omega_5_positive_20001bp.csv

cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/positive/fmh_omega_7.csv \
   figure2_datasets/fmh_omega_7_positive_5001bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/fmh_dnds_sketch_protein/positive_selection_redo_sketch_protein_using_faa/dnds.k7.approximations_included.csv \
   figure2_datasets/fmh_omega_7_positive_10002bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/positive/fmh_omega_7.csv \
   figure2_datasets/fmh_omega_7_positive_20001bp.csv

cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/positive/k15/fmh_omega_15.csv \
   figure2_datasets/fmh_omega_15_positive_5001bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/ksizes/15/positive/fmh_omega_15.csv \
   figure2_datasets/fmh_omega_15_positive_10002bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/positive/k15/fmh_omega_15.csv \
   figure2_datasets/fmh_omega_15_positive_20001bp.csv

cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/positive/k21/fmh_omega_21.csv \
   figure2_datasets/fmh_omega_21_positive_5001bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/ksizes/21/positive/fmh_omega_21.csv \
   figure2_datasets/fmh_omega_21_positive_10002bp.csv
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/positive/k21/fmh_omega_21.csv \
   figure2_datasets/fmh_omega_21_positive_20001bp.csv

# NG86 datasets (keep original names)
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/kaks_NG/positive_selection_queries_5001_0.01.axt.kaks \
   figure2_datasets/
cp /data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/kaks_NG/positive_selection_queries_10002_0.01.axt.kaks \
   figure2_datasets/
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/kaks_NG/positive_selection_queries_20001_0.01.axt.kaks \
   figure2_datasets/

cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/kaks_NG/negative_selection_queries_5001_0.01.axt.kaks \
   figure2_datasets/
cp /data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/kaks_NG/negative_selection_queries_10002_0.01.axt.kaks \
   figure2_datasets/
cp /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/kaks_NG/negative_selection_queries_20001_0.01.axt.kaks \
   figure2_datasets/

# Panel B (LAMA3) datasets
cp /data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/positive/fmh_k7/fmh_omega_7_modified.csv \
   figure2_datasets/fmh_omega_7_positive_LAMA3.csv

cp /data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/positive/fmh_k9/fmh_omega_9.csv \
   figure2_datasets/fmh_omega_9_positive_LAMA3.csv

cp /data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/positive/fmh_k11/fmh_omega_11.csv \
   figure2_datasets/fmh_omega_11_positive_LAMA3.csv

cp /data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/negative/fmh_k7/fmh_omega_7_modified.csv \
   figure2_datasets/fmh_omega_7_negative_LAMA3.csv

cp /data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/negative/fmh_k9/fmh_omega_9.csv \
   figure2_datasets/fmh_omega_9_negative_LAMA3.csv

cp /data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/negative/fmh_k11/fmh_omega_11.csv \
   figure2_datasets/fmh_omega_11_negative_LAMA3.csv

cp /data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/positive/kaks/kaks_sequences.axt.kaks \
   figure2_datasets/kaks_sequences_positive_LAMA3.axt.kaks

cp /data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/negative/kaks/kaks_sequences.axt.kaks \
   figure2_datasets/kaks_sequences_negative_LAMA3.axt.kaks