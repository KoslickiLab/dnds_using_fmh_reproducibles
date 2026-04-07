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

######### FILES FOR POSITIVE NG86 DNDS
ng86_dnds_5001_negative_files = '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/kaks_NG/negative_selection_queries_5001_0.01.axt.kaks'
ng86_dnds_10002_negative_files = '/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/kaks_NG/negative_selection_queries_10002_0.01.axt.kaks'
ng86_dnds_20001_negative_files = '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/kaks_NG/negative_selection_queries_20001_0.01.axt.kaks'

############ STORING FILES
fmh_dnds_negative_files = [fmh_dnds_k5_negative_files,
                    fmh_dnds_k7_negative_files,
                    fmh_dnds_k15_negative_files,
                    fmh_dnds_k21_negative_files,]

ng_dnds_negative_files = [ng86_dnds_5001_negative_files,
                            ng86_dnds_10002_negative_files,
                            ng86_dnds_20001_negative_files]

############ PLOTTING

# SUBPLOT
fig, axes = plt.subplots(3, 4,sharex=False, sharey=False, figsize=(64, 48))
fs=8.5
ticks_lst=[0,1]
limits=(-0.1, 1.5)

for i in range(3):
        axes[i,0].grid(False)
        axes[i,1].grid(False)
        axes[i,2].grid(False)
        axes[i,3].grid(False)

        axes[i,0].set_ylabel('NG86 dN/dS estimates',fontsize=fs)

        axes[i,1].set_yticks(ticks_lst)
        axes[i,0].set_yticks(ticks_lst)
        axes[i,2].set_yticks(ticks_lst)
        axes[i,3].set_yticks(ticks_lst)

        axes[i, 0].set_ylim(limits)
        axes[i, 1].set_ylim(limits)
        axes[i, 2].set_ylim(limits)
        axes[i, 3].set_ylim(limits)

        axes[i, 0].set_xlim(limits)
        axes[i, 1].set_xlim(limits)
        axes[i, 2].set_xlim(limits)
        axes[i, 3].set_xlim(limits)

        #RIGHT Y-AXES LABELS and TICKS
        axes[i,3].yaxis.set_label_position("right")

        axes[i,1].set_xticks(ticks_lst)
        axes[i,0].set_xticks(ticks_lst)
        axes[i,2].set_xticks(ticks_lst)
        axes[i,3].set_xticks(ticks_lst)

        axes[i,0].tick_params(axis='both', which='major', labelsize=fs-2)
        axes[i,1].tick_params(axis='both', which='major', labelsize=fs-2)
        axes[i,2].tick_params(axis='both', which='major', labelsize=fs-2)
        axes[i,3].tick_params(axis='both', which='major', labelsize=fs-2)

axes[2,0].set_xlabel('FMH dN/dS estimates',fontsize=fs)
axes[2,1].set_xlabel('FMH dN/dS estimates',fontsize=fs)
axes[2,2].set_xlabel('FMH dN/dS estimates',fontsize=fs)
axes[2,3].set_xlabel('FMH dN/dS estimates',fontsize=fs)

axes[0,0].set_title(r"$\it{k}$"+"=5, 5000 bp", fontsize=fs)
axes[0,1].set_title(r"$\it{k}$"+"=7, 5000 bp", fontsize=fs)
axes[0,2].set_title(r"$\it{k}$"+"=15, 5000 bp", fontsize=fs)
axes[0,3].set_title(r"$\it{k}$"+"=21, 5000 bp", fontsize=fs)

axes[1,0].set_title(r"$\it{k}$"+"=5, 10,002 bp", fontsize=fs)
axes[1,1].set_title(r"$\it{k}$"+"=7, 10,002 bp", fontsize=fs)
axes[1,2].set_title(r"$\it{k}$"+"=15, 10,002 bp", fontsize=fs)
axes[1,3].set_title(r"$\it{k}$"+"=21, 10,002 bp", fontsize=fs)

axes[2,0].set_title(r"$\it{k}$"+"=5, 20,001 bp", fontsize=fs)
axes[2,1].set_title(r"$\it{k}$"+"=7, 20,001 bp", fontsize=fs)
axes[2,2].set_title(r"$\it{k}$"+"=15, 20,001 bp", fontsize=fs)
axes[2,3].set_title(r"$\it{k}$"+"=21, 20,001 bp", fontsize=fs)

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

for length in range(3):
        for k in range(4):
                fmh_negative=pd.read_csv(fmh_dnds_negative_files[k][length],sep=',')
                fmh_negative=fmh_negative[fmh_negative['A']=='ref_gene'][['B','dN/dS']]
                fmh_negative=fmh_negative[(fmh_negative['dN/dS'] > 0) & (fmh_negative['dN/dS'] < 10)]

                ng86_negative=pd.read_csv(ng_dnds_negative_files[length],sep='\t')[['Sequence','Ka/Ks']]
                ng86_negative=ng86_negative[(ng86_negative['Ka/Ks'] > 0) & (ng86_negative['Ka/Ks'] < 10)]

                negative_df=pd.concat([fmh_negative,ng86_negative],axis=1).dropna()
                negative_corr, negative_corr_p = scipy.stats.pearsonr(negative_df['Ka/Ks'], negative_df['dN/dS'])

                axes[length, k].scatter(negative_df['dN/dS'], negative_df['Ka/Ks'], color='blue',label='Negative', alpha=0.5)
                axes[length, k].text(0.5, 1.25, f'Pearson R: {round(negative_corr,3)}',fontsize=fs-2)


plt.subplots_adjust(left=0.9,
                    bottom=0.9, 
                    right=1, 
                    top=1, 
                    wspace=0.1, 
                    hspace=0.2)  # Increase right margin

fig.figure.savefig(f"/data/jzr5814/repositories/dnds_using_fmh_reproducibles/manuscript_figures/updated_pdf/figure2a_sup_negative.pdf",bbox_inches='tight') 
