import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set font to Times New Roman
rcParams['font.family'] = 'serif'


######### FILES FOR NEGATIVE FMH DNDS
fmh_dnds_k9_negative_files = ['/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k9_10scale_cores1_includes_ref/fmh_omega_9.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k9_100scale_cores1_includes_ref/fmh_omega_9.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k9_1000scale_cores1_includes_ref/fmh_omega_9.csv']

######### FILES FOR POSITIVE FMH DNDS
fmh_dnds_k9_positive_files = ['/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k9_10scale_cores1_includes_ref/fmh_omega_9.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k9_100scale_cores1_includes_ref/fmh_omega_9.csv',
        '/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k9_1000scale_cores1_includes_ref/fmh_omega_9.csv']

nega1='/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k9_10scale_cores1_includes_ref/fmh_omega_9.csv'
nega1_df=pd.read_csv(nega1,sep=',')
nega1_df=nega1_df[nega1_df['query_name']=='ref_ecoli'][['match_name','dN/dS']]
nega1_df = nega1_df[(nega1_df['dN/dS'] >= 0) & (nega1_df['dN/dS'] <= 10)]


nega2='/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k9_100scale_cores1_includes_ref/fmh_omega_9.csv'
nega2_df=pd.read_csv(nega2,sep=',')
nega2_df=nega2_df[nega2_df['query_name']=='ref_ecoli'][['match_name','dN/dS']]
nega2_df = nega1_df[(nega2_df['dN/dS'] >= 0) & (nega2_df['dN/dS'] <= 10)]


nega3='/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k9_1000scale_cores1_includes_ref/fmh_omega_9.csv'
nega3_df=pd.read_csv(nega3,sep=',')
nega3_df=nega3_df[nega3_df['query_name']=='ref_ecoli'][['match_name','dN/dS']]
nega3_df = nega3_df[(nega1_df['dN/dS'] >= 0) & (nega3_df['dN/dS'] <= 10)]


posi1='/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k9_10scale_cores1_includes_ref/fmh_omega_9.csv'
posi1_df=pd.read_csv(posi1,sep=',')
posi1_df=posi1_df[posi1_df['query_name']=='ref_ecoli'][['match_name','dN/dS']]
posi1_df = posi1_df[(posi1_df['dN/dS'] >= 0) & (posi1_df['dN/dS'] <= 10)]


posi2='/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k9_100scale_cores1_includes_ref/fmh_omega_9.csv'
posi2_df=pd.read_csv(posi2,sep=',')
posi2_df=posi2_df[posi2_df['query_name']=='ref_ecoli'][['match_name','dN/dS']]
posi2_df = posi2_df[(posi2_df['dN/dS'] >= 0) & (posi2_df['dN/dS'] <= 10)]

posi3='/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k9_1000scale_cores1_includes_ref/fmh_omega_9.csv'
posi3_df=pd.read_csv(posi3,sep=',')
posi3_df=posi3_df[posi3_df['query_name']=='ref_ecoli'][['match_name','dN/dS']]
posi3_df = posi3_df[(posi3_df['dN/dS'] >= 0) & (posi3_df['dN/dS'] <= 10)]


############ PLOTTING



# SUBPLOT

fig, axes = plt.subplots(1, 3,sharex=True, sharey=True, figsize=(6, 2))
fs=10
bins = 50

range_min = min(nega1_df['dN/dS'].min(), nega2_df['dN/dS'].min(), nega3_df['dN/dS'].min(),
                posi1_df['dN/dS'].min(), posi2_df['dN/dS'].min(), posi3_df['dN/dS'].min())
range_max = max(nega1_df['dN/dS'].max(), nega2_df['dN/dS'].max(), nega3_df['dN/dS'].max(),
                posi1_df['dN/dS'].max(), posi2_df['dN/dS'].max(), posi3_df['dN/dS'].max())


axes[0].hist(nega1_df['dN/dS'], color='blue',bins=bins, range=(range_min, range_max), alpha=0.5, edgecolor='blue')
axes[0].hist(posi1_df['dN/dS'], color='darkorange',bins=bins, range=(range_min, range_max), alpha=0.5, edgecolor='darkorange')
axes[0].tick_params(axis='both', which='major', labelsize=fs-2)

axes[1].hist(nega2_df['dN/dS'], color='blue',bins=bins, range=(range_min, range_max), alpha=0.5, edgecolor='blue')
axes[1].hist(posi2_df['dN/dS'], color='darkorange',bins=bins, range=(range_min, range_max), alpha=0.5, edgecolor='darkorange')
axes[1].tick_params(axis='both', which='major', labelsize=fs-2)

axes[2].hist(nega3_df['dN/dS'], color='blue',bins=bins, range=(range_min, range_max), alpha=0.5, edgecolor='blue')
axes[2].hist(posi3_df['dN/dS'], color='darkorange',bins=bins, range=(range_min, range_max), alpha=0.5, edgecolor='darkorange')
axes[2].tick_params(axis='both', which='major', labelsize=fs-2)

scales=[10,100,1000]
for scalef in range(3):
        axes[scalef].axvline(x=1, linestyle='--', color='grey')
        axes[scalef].set_xlabel('FMH dN/dS estimates',fontsize=fs)
        axes[scalef].set_title(f"scaled={scales[scalef]}",fontsize=fs)
        axes[scalef].set_xlim(-0.5,8.5)
        axes[scalef].set_xticks([0, 2, 4, 6, 8])

# Adjust layout
plt.subplots_adjust(right=1)  # Increase right margin
plt.subplots_adjust(wspace=0.1, hspace=0.15)
fig.figure.savefig(f"/data/jzr5814/sourmash_dnds_estimation/thesis_figures/sup2.pdf",bbox_inches='tight') 
