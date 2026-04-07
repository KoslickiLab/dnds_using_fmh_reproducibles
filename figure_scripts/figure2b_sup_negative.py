import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
from matplotlib import rcParams

# Set font to Times New Roman
rcParams['font.family'] = 'serif'

######### FILES FOR NEGATIVE FMH DNDS
fmh_dnds_negative_files = ['/data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/negative/fmh_k7/fmh_omega_7_modified.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/negative/fmh_k9/fmh_omega_9.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/negative/fmh_k11/fmh_omega_11.csv']

######### FILES FOR POSITIVE NG86 DNDS
ng86_dnds_negative_files = '/data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/negative/kaks/kaks_sequences.axt.kaks'

############ PLOTTING

# SUBPLOT
fig, axes = plt.subplots(1, 3,sharex=True, sharey=True, figsize=(6, 2))
fs=8.5
fs_label=fs
ticks_lst=[0,1]
limits=(-0.1, 1.75)

# AXES
for i in range(3):
        axes[i].sharex(axes[0])
        axes[i].sharey(axes[0])

for i in range(3):

        axes[i].set_yticks(ticks_lst)
        axes[i].set_ylim(limits)
        axes[i].set_xlim(limits)  
        axes[i].set_xticks(ticks_lst)
        axes[i].tick_params(axis='both', which='major', labelsize=fs-2)
        axes[i].set_xlabel('FMH dN/dS estimates',fontsize=fs)

axes[0].set_title(r"$\it{k}$"+"=7",fontsize=fs)
axes[1].set_title(r"$\it{k}$"+"=9",fontsize=fs)
axes[2].set_title(r"$\it{k}$"+"=11",fontsize=fs)

axes[0].set_ylabel('NG86 dN/dS estimates',fontsize=fs)

for k in range(3):
        fmh_negative=pd.read_csv(fmh_dnds_negative_files[k],sep=',')
        fmh_negative=fmh_negative[fmh_negative['A']=='negative_0.01_ref'][['B','dN/dS']]
        fmh_negative=fmh_negative[(fmh_negative['dN/dS'] > 0) & (fmh_negative['dN/dS'] < 8)]

        ng86_negative=pd.read_csv(ng86_dnds_negative_files,sep='\t')[['Sequence','Ka/Ks']]
        ng86_negative=ng86_negative[ng86_negative['Sequence']=='negative_0.01_ref']
        ng86_negative=ng86_negative[(ng86_negative['Ka/Ks'] > 0) & (ng86_negative['Ka/Ks'] < 8)]

        negative_df=pd.concat([fmh_negative,ng86_negative],axis=1).dropna()

        negative_corr, negative_corr_p = scipy.stats.pearsonr(negative_df['Ka/Ks'], negative_df['dN/dS'])


        axes[k].scatter(negative_df['dN/dS'], negative_df['Ka/Ks'], color='blue',label='Negative', alpha=0.5)
        axes[k].text(0.75, 1.5, f'Pearson R: {round(negative_corr,3)}',fontsize=fs-2)
        

plt.subplots_adjust(right=1)  # Increase right margin
plt.subplots_adjust(wspace=0.1, hspace=0.15)
fig.figure.savefig(f"/data/jzr5814/repositories/dnds_using_fmh_reproducibles/manuscript_figures/updated_pdf/figure2b_sup_negative.pdf",bbox_inches='tight') 
