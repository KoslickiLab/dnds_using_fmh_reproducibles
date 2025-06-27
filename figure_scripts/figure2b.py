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


######### FILES FOR POSITIVE FMH DNDS
fmh_dnds_positive_files = ['/data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/positive/fmh_k7/fmh_omega_7_modified.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/positive/fmh_k9/fmh_omega_9.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/positive/fmh_k11/fmh_omega_11.csv']

######### FILES FOR POSITIVE NG86 DNDS
ng86_dnds_positive_files = '/data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/positive/kaks/kaks_sequences.axt.kaks'

######### FILES FOR POSITIVE NG86 DNDS
ng86_dnds_negative_files = '/data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/negative/kaks/kaks_sequences.axt.kaks'

############ PLOTTING

# SUBPLOT
fig, axes = plt.subplots(1, 3,sharex=True, sharey=True, figsize=(6, 2))
fs=10
fs_label=fs

# AXES
for i in range(3):
        axes[i].sharex(axes[0])
        axes[i].sharey(axes[0])

for i in range(3):

        axes[i].set_yticks([0,1,2,5,8])
        axes[i].set_ylim(0,8)
        axes[i].set_xlim(0,8)  
        axes[i].set_xticks([0,1,2,5,8])
        axes[i].tick_params(axis='both', which='major', labelsize=fs-2)
        axes[i].set_xlabel('FMH dN/dS estimates',fontsize=fs)

axes[0].set_title(r"$\it{k}$"+"=7",fontsize=fs)
axes[1].set_title(r"$\it{k}$"+"=9",fontsize=fs)
axes[2].set_title(r"$\it{k}$"+"=11",fontsize=fs)

axes[0].set_ylabel('NG86 dN/dS estimates',fontsize=fs)

for k in range(3):
        fmh_positive=pd.read_csv(fmh_dnds_positive_files[k],sep=',')
        fmh_positive=fmh_positive[fmh_positive['A']=='positive_0.01_ref'][['B','dN/dS']]
        fmh_positive=fmh_positive[(fmh_positive['dN/dS'] > 0) & (fmh_positive['dN/dS'] < 8)]
        #print('fmh_positive',fmh_positive.describe())
        
        fmh_negative=pd.read_csv(fmh_dnds_negative_files[k],sep=',')
        fmh_negative=fmh_negative[fmh_negative['A']=='negative_0.01_ref'][['B','dN/dS']]
        fmh_negative=fmh_negative[(fmh_negative['dN/dS'] > 0) & (fmh_negative['dN/dS'] < 8)]
        #print('fmh_negative',fmh_negative.describe())

        ng86_positive=pd.read_csv(ng86_dnds_positive_files,sep='\t')[['Sequence','Ka/Ks']]
        ng86_positive=ng86_positive[ng86_positive['Sequence']=='positive_0.01_ref']
        ng86_positive=ng86_positive[(ng86_positive['Ka/Ks'] > 0) & (ng86_positive['Ka/Ks'] < 8)]


        ng86_negative=pd.read_csv(ng86_dnds_negative_files,sep='\t')[['Sequence','Ka/Ks']]
        ng86_negative=ng86_negative[ng86_negative['Sequence']=='negative_0.01_ref']
        ng86_negative=ng86_negative[(ng86_negative['Ka/Ks'] > 0) & (ng86_negative['Ka/Ks'] < 8)]
        #print(ng86_negative.describe())

        positive_df=pd.concat([fmh_positive,ng86_positive],axis=1).dropna()
        negative_df=pd.concat([fmh_negative,ng86_negative],axis=1).dropna()

        combined = pd.concat([positive_df, negative_df], ignore_index=True)
        combined_corr, combined_corr_p = scipy.stats.pearsonr(combined['Ka/Ks'], combined['dN/dS'])

        axes[k].scatter(negative_df['dN/dS'], negative_df['Ka/Ks'], color='blue',label='Negative', alpha=0.5)
        axes[k].scatter(positive_df['dN/dS'], positive_df['Ka/Ks'], color='darkorange',label='Positive', alpha=0.5)

        tmp = [negative_df['Ka/Ks'].min(), positive_df['Ka/Ks'].max()]
        axes[k].plot([tmp[0], tmp[1]], [tmp[0], tmp[1]], linestyle='--', color='grey')
        axes[k].text(3, 7.3, f'Pearson R: {round(combined_corr,3)}',fontsize=fs-2)
        

plt.subplots_adjust(right=1)  # Increase right margin
plt.subplots_adjust(wspace=0.1, hspace=0.15)
fig.figure.savefig(f"/data/jzr5814/sourmash_dnds_estimation/thesis_figures/figure2b.pdf",bbox_inches='tight') 
