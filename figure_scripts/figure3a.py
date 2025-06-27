import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
from matplotlib import rcParams

# Set font to Times New Roman
rcParams['font.family'] = 'serif'

#create figure that compares interpretations to martinez-gutierrez(2022)
ksize=7
scaled=500

PdS_threshold=0.05
PdN_threshold=0.23

thresh=0.02
data1=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_{scaled}scale_{thresh}threshold_k{ksize}/fmh_omega_{ksize}_by_genome_group_dSminthresh{PdS_threshold}_and_dNminthresh{PdN_threshold}.csv'
data_k7_t002=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k7/fmh_omega_7_by_genome_group_dSminthresh0.05_and_dNminthresh0.23.csv'

thresh=0.05
data_k7_t005=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k7/fmh_omega_7_by_genome_group_dSminthresh0.05_and_dNminthresh0.23.csv'

thresh=0.1
data_k7_t010=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k7/fmh_omega_7_by_genome_group_dSminthresh0.05_and_dNminthresh0.23.csv'

thresh=0.2
data_k7_t020=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_500scale_{thresh}threshold_k7/fmh_omega_7_by_genome_group_dSminthresh0.05_and_dNminthresh0.23.csv'

paper='/data/jzr5814/sourmash_dnds_estimation/tests/data/GTDB/martinez-gutierrez_2022paper/pgen.1010220.s002.tsv'
paper_df=pd.read_csv(f'{paper}',sep='\t')[['Genome','Genus','Mediandnds','Meandnds']].set_index('Genome')
tmp = [min(paper_df['Mediandnds']), max(paper_df['Mediandnds'])]

data={"data_k7_t0.02":{"filename":data_k7_t002},
      "data_k7_t0.05":{"filename":data_k7_t005},
      "data_k7_t0.1":{"filename":data_k7_t010},
      "data_k7_t0.2":{"filename":data_k7_t020}}

for file in data:
    df = pd.read_csv(data[file]["filename"],sep=',')[['mean', 'median', 'genome']]
    df['genome'] = df['genome'].apply(lambda row: '_'.join(row.split('_')[1:]) if '_' in row else row)
    df = df.set_index('genome')
    #df=df[df["median"]<=0.4]

    results = pd.concat([df,paper_df],axis=1).dropna()
    median_corr, median_ =scipy.stats.pearsonr(results['median'], results['Mediandnds'])

    data[file]["results"] = results
    data[file]["median_dnds_r"] = median_corr

####### Plot

fig, axes = plt.subplots(1, 4,sharey=True,sharex=True,figsize=(16,4))
#plt.subplots_adjust(wspace=0.01,
#                    hspace=0.15)

fs=15
thresholds = [0.02, 0.05, 0.1, 0.2]

for t in range(len(thresholds)):
    axes[t].scatter(data[f"data_k7_t{thresholds[t]}"]["results"]['median'], data[f"data_k7_t{thresholds[t]}"]["results"]['Mediandnds'], color='blue', 
        edgecolor='blue', alpha=0.5)
    axes[t].set_xlabel('FMH dN/dS estimates',fontsize=fs)

for t in range(len(thresholds)):
    for x in range(3):
        axes[t].plot(tmp, tmp, linestyle='--', color='black', alpha=1)
        axes[t].set_xlim(-0.04,0.5)
        axes[t].set_ylim(-0.04,0.5)
        axes[t].set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
        axes[t].set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
        axes[t].tick_params(axis='x', labelsize=fs-2)  
        axes[t].tick_params(axis='y', labelsize=fs-2)  
        axes[t].set_title(f"t={thresholds[t]},k=7",fontsize=fs)
        correlation = data[f"data_k7_t{thresholds[t]}"]['median_dnds_r']
        axes[t].text(0.2, 0.4, f'Pearson R: {round(correlation,3)}', fontsize = fs-2)
        axes[t].text(0.21, 0.21, f'x=y', fontsize = fs-2)

axes[0].set_ylabel('CodeML dN/dS estimates',fontsize=fs)

# Adjust layout
plt.tight_layout()
fig.figure.savefig(f"/data/jzr5814/sourmash_dnds_estimation/thesis_figures/figure3a.pdf",bbox_inches='tight') 
