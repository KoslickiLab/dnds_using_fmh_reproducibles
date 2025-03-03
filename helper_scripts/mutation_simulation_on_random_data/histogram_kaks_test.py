#!/usr/bin/env python3

'''
    This script creates a figure to evaluate the selection simulation on random data 
    to check whether KaKs_calculator using the NG86 model can indicate positive and negative selection. 
    If there is a distinction between both selections, then that indicates the simulation works 
    and move on analyzing FMH dn/ds estimations on random data.
'''

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set font to Times New Roman
rcParams['font.family'] = 'serif'

# Create subplots
fig, axes = plt.subplots(1, 3,sharex=True, sharey=True, figsize=(14, 5))
plt.subplots_adjust(right=1,
                    wspace=0.08, 
                    hspace=0.15)  # Increase right margin
fs=15
bin_size=50
xlabel="NG86 dN/dS estimates"

### plotting analysis for p_rate = 0.1

wd='/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.1/kaks_NG'
file1=pd.read_csv(f'{wd}/negative_selection_queries_10002_0.1.axt.kaks',sep='\t').set_index("Sequence")['Ka/Ks']
file2=pd.read_csv(f'{wd}/positive_selection_queries_10002_0.1.axt.kaks',sep='\t').set_index("Sequence")['Ka/Ks']
#posianifile = pd.read_csv('/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.1/fmh_dnds/positive_selection/ani_k21_against_ref_gene.csv').set_index("Unnamed: 0")[["0"]]
#negaanifile = pd.read_csv('/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.1/fmh_dnds/negative_selection/ani_k21_against_ref_gene.csv').set_index("Unnamed: 0")[["0"]]
#result_negative = pd.concat([file1, negaanifile], axis=1)
result_negative = result_negative[(result_negative['Ka/Ks'] >= 0) & (result_negative['Ka/Ks'] <= 10)]

result_positive = pd.concat([file2, posianifile], axis=1)
result_positive = result_positive[(result_positive['Ka/Ks'] >= 0) & (result_positive['Ka/Ks'] <= 10)]

# Plot the first histogram with primary y-axis
axes[0].hist(result_negative['Ka/Ks'], edgecolor='blue',alpha=0.5, label='negative',bins=bin_size)
axes[0].hist(result_positive['Ka/Ks'], edgecolor='darkorange',alpha=0.5, label='positive',bins=bin_size)
axes[0].axvline(x=1, color='gray', linestyle='--', linewidth=1)
axes[0].set_xlabel(f'{xlabel}',fontsize=fs-2)
axes[0].set_title(r'$\it{p}$'+'=0.1',fontsize=fs)
axes[0].grid(False)
axes[0].tick_params(axis='both', labelsize=fs-3)

### plotting analysis for p_rate = 0.01

wd='/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/kaks_NG'
file1=pd.read_csv(f'{wd}/negative_selection_queries_10002_0.01.axt.kaks',sep='\t').set_index("Sequence")['Ka/Ks']
file2=pd.read_csv(f'{wd}/positive_selection_queries_10002_0.01.axt.kaks',sep='\t').set_index("Sequence")['Ka/Ks']
#posianifile = pd.read_csv('/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/fmh_dnds_sketch_translate/positive_selection/ani_k21_against_ref_gene.csv').set_index("Unnamed: 0")[["0"]]
#negaanifile = pd.read_csv('/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/fmh_dnds_sketch_translate/negative_selection/ani_k21_against_ref_gene.csv').set_index("Unnamed: 0")[["0"]]
#result_negative = pd.concat([file1, negaanifile], axis=1)
result_negative = result_negative[(result_negative['Ka/Ks'] >= 0) & (result_negative['Ka/Ks'] <= 10)]

result_positive = pd.concat([file2, posianifile], axis=1)
result_positive = result_positive[(result_positive['Ka/Ks'] >= 0) & (result_positive['Ka/Ks'] <= 10)]

axes[1].hist(result_negative['Ka/Ks'], edgecolor='blue', alpha=0.5, label='negative',bins=bin_size)
axes[1].hist(result_positive['Ka/Ks'], edgecolor='darkorange', alpha=0.5, label='positive',bins=bin_size)
axes[1].axvline(x=1, color='gray', linestyle='--', linewidth=1)
axes[1].set_xlabel(f'{xlabel}',fontsize=fs-2)
axes[1].set_title(r'$\it{p}$'+'=0.01',fontsize=fs)
axes[1].grid(False)
axes[1].tick_params(axis='both', labelsize=fs-3)

### plotting analysis for p_rate = 0.001

wd='/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.001/kaks_NG'
file1=pd.read_csv(f'{wd}/negative_selection_queries_10002_0.001.axt.kaks',sep='\t').set_index("Sequence")['Ka/Ks']
file2=pd.read_csv(f'{wd}/positive_selection_queries_10002_0.001.axt.kaks',sep='\t').set_index("Sequence")['Ka/Ks']
#posianifile = pd.read_csv('/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.001/fmh_dnds_sketch_translate/positive_selection/ani_k21_against_ref_gene.csv').set_index("Unnamed: 0")[["0"]]
#negaanifile = pd.read_csv('/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.001/fmh_dnds_sketch_translate/negative_selection/ani_k21_against_ref_gene.csv').set_index("Unnamed: 0")[["0"]]
#result_negative = pd.concat([file1, negaanifile], axis=1)
result_negative = result_negative[(result_negative['Ka/Ks'] >= 0) & (result_negative['Ka/Ks'] <= 10)]

result_positive = pd.concat([file2, posianifile], axis=1)
result_positive = result_positive[(result_positive['Ka/Ks'] >= 0) & (result_positive['Ka/Ks'] <= 10)]

axes[2].hist(result_negative['Ka/Ks'], edgecolor='blue', alpha=0.5, label='negative',bins=bin_size)
axes[2].hist(result_positive['Ka/Ks'], edgecolor='darkorange', alpha=0.5, label='positive',bins=bin_size)
axes[2].axvline(x=1, color='gray', linestyle='--', linewidth=1)
axes[2].set_xlabel(f'{xlabel}',fontsize=fs-2)
axes[2].set_title(r'$\it{p}$'+'=0.001',fontsize=fs)
axes[2].grid(False)
axes[2].tick_params(axis='both', labelsize=fs-3)

### Save figure as pdf

fig.figure.savefig(f"/data/jzr5814/sourmash_dnds_estimation/thesis_figures/histogram_kaks_test_v7.pdf",bbox_inches='tight') 


