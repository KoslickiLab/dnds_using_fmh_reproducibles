#!/usr/bin/env python3

'''
    This script creates figure of pairwise pearson correlation fmh dnds estimations across different parameters 
    (i.e., ksizes and containment thresholds)
'''

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
from matplotlib import rcParams
import seaborn as sns
 
###### PREPROCESS DATA

#create figure that represents median Cfracs across ksizes for both DNA and Protein
thresh=0.02
data_k7_t002=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k7/fmh_omega_7.csv'
data_k9_t002=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k9/fmh_omega_9.csv'
data_k11_t002=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k11/fmh_omega_11.csv'

thresh=0.05
data_k7_t005=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k7/fmh_omega_7.csv'
data_k9_t005=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k9/fmh_omega_9.csv'
data_k11_t005=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k11/fmh_omega_11.csv'

thresh=0.1
data_k7_t010=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k7/fmh_omega_7.csv'
data_k9_t010=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k9/fmh_omega_9.csv'
data_k11_t010=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k11/fmh_omega_11.csv'

thresh=0.2
data_k7_t020=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k7/fmh_omega_7.csv'
data_k9_t020=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k9/fmh_omega_9.csv'
data_k11_t020=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k11/fmh_omega_11.csv'

data={"k_aa=7, t=0.02":{"filename":data_k7_t002}, "k_aa=9, t=0.02":{"filename":data_k9_t002}, "k_aa=11, t=0.02":{"filename":data_k11_t002},
      "k_aa=7, t=0.05":{"filename":data_k7_t005}, "k_aa=9, t=0.05":{"filename":data_k9_t005}, "k_aa=11, t=0.05":{"filename":data_k11_t005},
      "k_aa=7, t=0.1":{"filename":data_k7_t010}, "k_aa=9, t=0.1":{"filename":data_k9_t010}, "k_aa=11, t=0.1":{"filename":data_k11_t010},
      "k_aa=7, t=0.2":{"filename":data_k7_t020}, "k_aa=9, t=0.2":{"filename":data_k9_t020}, "k_aa=11, t=0.2":{"filename":data_k11_t020}}

df_list = []

for file in data:
    temp_df = pd.read_csv(data[file]["filename"],sep=',')[["A,B",'dN/dS']].set_index("A,B").rename(columns={"dN/dS":file})
    temp_df=temp_df[(temp_df[file]>0) & (temp_df[file]<=2)]
    df_list.append(temp_df)

df = pd.concat(df_list, axis=1).dropna()

# Create the pairplot
fig = sns.pairplot(df,hue=None,
                    plot_kws={'color': 'gray', 'edgecolor': 'black'},
                    diag_kws={'color': 'gray', 'edgecolor': 'black', 'bins': 100})

# Apply the same xlim and ylim to all plots
for ax in fig.axes.flatten():
    if ax is not None:  # Some axes might be empty in pairplot
        ax.set_xlim(0, 2)
        ax.set_ylim(0, 2)

# Save the figure
fig.savefig("pairplot.png")  