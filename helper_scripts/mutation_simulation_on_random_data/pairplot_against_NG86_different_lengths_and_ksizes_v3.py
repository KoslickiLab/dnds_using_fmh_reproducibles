#!/usr/bin/env python3

'''
    This script creates a figure to evaluate FMH dn/ds estimations on random data 
    to check whether FMH dn/ds estimations can distinguish between negative and positive selections across different ksizes and sequence lengths
    and move on analyzing FMH dn/ds estimations on simulated mutation on real data.
'''

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
from matplotlib import rcParams

# Set font to Times New Roman
rcParams['font.family'] = 'serif'

######### FILES FOR NEGATIVE FMH DNDS
fmh_dnds_k5_negative_files = ['/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/negative/k5/fmh_omega_5.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/ksizes/5/negative/fmh_omega_5.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/negative/k5/fmh_omega_5.csv']

fmh_dnds_k7_negative_files = ['/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/negative/fmh_omega_7.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/fmh_dnds_sketch_protein/negative_selection_redo_sketch_protein_using_faa/dnds.k7.approximations_included.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/negative/fmh_omega_7.csv']

fmh_dnds_k15_negative_files = ['/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/negative/k15/fmh_omega_15.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/ksizes/15/negative/fmh_omega_15.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/negative/k15/fmh_omega_15.csv']

fmh_dnds_k21_negative_files = ['/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/negative/k21/fmh_omega_21.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/ksizes/21/negative/fmh_omega_21.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/negative/k21/fmh_omega_21.csv']


######### FILES FOR POSITIVE FMH DNDS
fmh_dnds_k5_positive_files = ['/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/positive/k5/fmh_omega_5.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/ksizes/5/positive/fmh_omega_5.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/positive/k5/fmh_omega_5.csv']

fmh_dnds_k7_positive_files = ['/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/positive/fmh_omega_7.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/fmh_dnds_sketch_protein/positive_selection_redo_sketch_protein_using_faa/dnds.k7.approximations_included.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/positive/fmh_omega_7.csv']

fmh_dnds_k15_positive_files = ['/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/positive/k15/fmh_omega_15.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/ksizes/15/positive/fmh_omega_15.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/positive/k15/fmh_omega_15.csv']

fmh_dnds_k21_positive_files = ['/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/positive/k21/fmh_omega_21.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/ksizes/21/positive/fmh_omega_21.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/positive/k21/fmh_omega_21.csv']


######### FILES FOR POSITIVE NG86 DNDS
ng86_dnds_5001_positive_files = '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/kaks_NG/positive_selection_queries_5001_0.01.axt.kaks'
ng86_dnds_10002_positive_files = '/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/kaks_NG/positive_selection_queries_10002_0.01.axt.kaks'
ng86_dnds_20001_positive_files = '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/kaks_NG/positive_selection_queries_20001_0.01.axt.kaks'

######### FILES FOR POSITIVE NG86 DNDS
ng86_dnds_5001_negative_files = '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/kaks_NG/negative_selection_queries_5001_0.01.axt.kaks'
ng86_dnds_10002_negative_files = '/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/kaks_NG/negative_selection_queries_10002_0.01.axt.kaks'
ng86_dnds_20001_negative_files = '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/kaks_NG/negative_selection_queries_20001_0.01.axt.kaks'

############ STORING FILES
fmh_dnds_negative_files = [fmh_dnds_k5_negative_files,
                    fmh_dnds_k7_negative_files,
                    fmh_dnds_k15_negative_files,
                    fmh_dnds_k21_negative_files,]

fmh_dnds_positive_files = [fmh_dnds_k5_positive_files,
                    fmh_dnds_k7_positive_files,
                    fmh_dnds_k15_positive_files,
                    fmh_dnds_k21_positive_files,]

ng_dnds_negative_files = [ng86_dnds_5001_negative_files,
                            ng86_dnds_10002_negative_files,
                            ng86_dnds_20001_negative_files]

ng_dnds_positive_files = [ng86_dnds_5001_positive_files,
                            ng86_dnds_10002_positive_files,
                            ng86_dnds_20001_positive_files]



############ PLOTTING

# SUBPLOT
scale_factor = 5
fig, axes = plt.subplots(3, 4,sharex=False, sharey=False, figsize=(50, 40))
fs=8.5

for i in range(3):
        axes[i,0].grid(False)
        axes[i,1].grid(False)
        axes[i,2].grid(False)
        axes[i,3].grid(False)

        axes[i,0].set_ylabel('NG86 dN/dS estimates',fontsize=fs-3)

        axes[i,1].set_yticks([0,1,2,5,10])
        axes[i,0].set_yticks([0,1,2,5,10])
        axes[i,2].set_yticks([0,1,2,5,10])
        axes[i,3].set_yticks([0,1,2,5,10])


        axes[i, 0].set_ylim(-1, 11)
        axes[i, 1].set_ylim(-1, 11)
        axes[i, 2].set_ylim(-1, 11)
        axes[i, 3].set_ylim(-1, 11)

        axes[i, 0].set_xlim(-1, 11)
        axes[i, 1].set_xlim(-1, 11)
        axes[i, 2].set_xlim(-1, 11)
        axes[i, 3].set_xlim(-1, 11)


        #RIGHT Y-AXES LABELS and TICKS
        axes[i,3].yaxis.set_label_position("right")

        axes[i,1].set_xticks([0,1,2,5,10])
        axes[i,0].set_xticks([0,1,2,5,10])
        axes[i,2].set_xticks([0,1,2,5,10])
        axes[i,3].set_xticks([0,1,2,5,10])

        axes[i,0].tick_params(axis='both', which='major', labelsize=fs-4)
        axes[i,1].tick_params(axis='both', which='major', labelsize=fs-4)
        axes[i,2].tick_params(axis='both', which='major', labelsize=fs-4)
        axes[i,3].tick_params(axis='both', which='major', labelsize=fs-4)



axes[2,0].set_xlabel('FMH Omega dN/dS estimates',fontsize=fs-3)
axes[2,1].set_xlabel('FMH Omega dN/dS estimates',fontsize=fs-3)
axes[2,2].set_xlabel('FMH Omega dN/dS estimates',fontsize=fs-3)
axes[2,3].set_xlabel('FMH Omega dN/dS estimates',fontsize=fs-3)


axes[0,0].set_title(r"$\it{k}$"+"=5", fontsize=fs-3)
axes[0,1].set_title(r"$\it{k}$"+"=7", fontsize=fs-3)
axes[0,2].set_title(r"$\it{k}$"+"=15", fontsize=fs-3)
axes[0,3].set_title(r"$\it{k}$"+"=21", fontsize=fs-3)

axes[0,0].set_xticklabels([])
axes[1,0].set_xticklabels([])

axes[0,1].set_xticklabels([])
axes[0,2].set_xticklabels([])
axes[0,3].set_xticklabels([])

axes[0,1].set_yticklabels([])
axes[0,2].set_yticklabels([])
axes[0,3].set_yticklabels([])


axes[1,1].set_xticklabels([])
axes[1,2].set_xticklabels([])
axes[1,3].set_xticklabels([])

axes[1,1].set_yticklabels([])
axes[1,2].set_yticklabels([])
axes[1,3].set_yticklabels([])

axes[2,1].set_yticklabels([])
axes[2,2].set_yticklabels([])
axes[2,3].set_yticklabels([])


rotate='horizontal'
#axes[0,3].set_ylabel("5,001 nt\n"+r"$\it{p}$=0.01",rotation=f'{rotate}', labelpad=35,fontsize=fs)
#axes[1,3].set_ylabel("10,002 nt\n"+r"$\it{p}$=0.01",rotation=f'{rotate}', labelpad=35,fontsize=fs)
#axes[2,3].set_ylabel("20,001 nt\n"+r"$\it{p}$=0.01",rotation=f'{rotate}', labelpad=35,fontsize=fs)

axes[0,3].set_ylabel("5,001 nt",rotation=f'{rotate}', labelpad=20,fontsize=fs-3)
axes[1,3].set_ylabel("10,002 nt",rotation=f'{rotate}', labelpad=20,fontsize=fs-3)
axes[2,3].set_ylabel("20,001 nt",rotation=f'{rotate}', labelpad=20,fontsize=fs-3)

#axes[0,0].set_ylabel('5,000 bp\nNG86 dN/dS estimates',fontsize=fs-3)
#axes[1,0].set_ylabel('10,002 bp\nNG86 dN/dS estimates',fontsize=fs-3)
#axes[2,0].set_ylabel('20,001 bp\nNG86 dN/dS estimates',fontsize=fs-3)


for length in range(3):
        for k in range(4):
                fmh_positive=pd.read_csv(fmh_dnds_positive_files[k][length],sep=',')
                fmh_positive=fmh_positive[fmh_positive['A']=='ref_gene'][['B','dN/dS']]
                fmh_positive=fmh_positive[(fmh_positive['dN/dS'] > 0) & (fmh_positive['dN/dS'] < 10)]

                fmh_negative=pd.read_csv(fmh_dnds_negative_files[k][length],sep=',')
                fmh_negative=fmh_negative[fmh_negative['A']=='ref_gene'][['B','dN/dS']]
                fmh_negative=fmh_negative[(fmh_negative['dN/dS'] > 0) & (fmh_negative['dN/dS'] < 10)]

                ng86_positive=pd.read_csv(ng_dnds_positive_files[length],sep='\t')[['Sequence','Ka/Ks']]
                ng86_positive=ng86_positive[(ng86_positive['Ka/Ks'] > 0) & (ng86_positive['Ka/Ks'] < 10)]

                ng86_negative=pd.read_csv(ng_dnds_negative_files[length],sep='\t')[['Sequence','Ka/Ks']]
                ng86_negative=ng86_negative[(ng86_negative['Ka/Ks'] > 0) & (ng86_negative['Ka/Ks'] < 10)]

                positive_df=pd.concat([fmh_positive,ng86_positive],axis=1).dropna()
                negative_df=pd.concat([fmh_negative,ng86_negative],axis=1).dropna()
                #nega_corr, nega_ = scipy.stats.pearsonr(negative_df['Ka/Ks'], negative_df['dN/dS'])
                #posi_corr, posi_ = scipy.stats.pearsonr(positive_df['Ka/Ks'], positive_df['dN/dS'])

                combined = pd.concat([positive_df, negative_df], ignore_index=True)
                combined_corr,combined_corr_p = scipy.stats.pearsonr(combined['Ka/Ks'], combined['dN/dS'])

                axes[length, k].scatter(negative_df['dN/dS'], negative_df['Ka/Ks'], color='blue',label='Negative', alpha=0.5)
                axes[length, k].scatter(positive_df['dN/dS'], positive_df['Ka/Ks'], color="darkorange",label='Positive', alpha=0.5)

                tmp = [negative_df['Ka/Ks'].min(), positive_df['Ka/Ks'].max()]

                axes[length, k].plot([tmp[0], tmp[1]], [tmp[0], tmp[1]], linestyle='--', color='grey')
                axes[length, k].text(5, 10, f'Pearson R: {"{:.3f}".format(combined_corr)}',fontsize=fs-4)
                #axes[k,length].text(1.75, 6.5, f'p-value: {round(posi_,3)}',fontsize=15)


# Adjust layout
#plt.suptitle(f'FMH dN/dS vs. NG86 dN/dS, k=7, p-rate=0.01, s=1',fontsize=15,y=0.95)
#fig.legend( loc='upper center', bbox_to_anchor=(0.5, 0.02), fancybox=True, shadow=True, ncol=2)
plt.subplots_adjust(left=0.9,
                    bottom=0.9, 
                    right=1, 
                    top=1, 
                    wspace=0.1, 
                    hspace=0.1)  # Increase right margin
fig.figure.savefig(f"/data/jzr5814/sourmash_dnds_estimation/thesis_figures/pairplot_against_NG86_scatter_lengths_and_ksizes_v12.pdf",bbox_inches='tight') 
