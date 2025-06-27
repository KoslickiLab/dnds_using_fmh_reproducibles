import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
from matplotlib import rcParams

# Set font to Times New Roman
rcParams['font.family'] = 'serif'


#create figure that represents median Cfracs across ksizes for both DNA and Protein
ksize=7
scaled=500
dSminthresh=0.05
dNminthresh=0.23

thresh=0.02
data1=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_{scaled}scale_{thresh}threshold_k{ksize}/fmh_omega_{ksize}_by_genus_group_dSminthresh{dSminthresh}_and_dNminthresh{dNminthresh}.csv'
data_df1=pd.read_csv(f'{data1}',sep=',')[['median','genus']].rename(columns={'median':'median1','genus':'Genus'}).set_index('Genus')

thresh=0.05
data2=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_{scaled}scale_{thresh}threshold_k{ksize}/fmh_omega_{ksize}_by_genus_group_dSminthresh{dSminthresh}_and_dNminthresh{dNminthresh}.csv'
data_df2=pd.read_csv(f'{data2}',sep=',')[['median','genus']].rename(columns={'median':'median2','genus':'Genus'}).set_index('Genus')

thresh=0.1
data3=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_{scaled}scale_{thresh}threshold_k{ksize}/fmh_omega_{ksize}_by_genus_group_dSminthresh{dSminthresh}_and_dNminthresh{dNminthresh}.csv'
data_df3=pd.read_csv(f'{data3}',sep=',')[['median','genus']].rename(columns={'median':'median3','genus':'Genus'}).set_index('Genus')

thresh=0.2
data4=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_{scaled}scale_{thresh}threshold_k{ksize}/fmh_omega_{ksize}_by_genus_group_dSminthresh{dSminthresh}_and_dNminthresh{dNminthresh}.csv'
data_df4=pd.read_csv(f'{data4}',sep=',')[['median','genus']].rename(columns={'median':'median3','genus':'Genus'}).set_index('Genus')

paper='/data/jzr5814/sourmash_dnds_estimation/tests/data/GTDB/martinez-gutierrez_2022paper/pgen.1010220.s002.tsv'
paper_df=pd.read_csv(f'{paper}',sep='\t')[['Genus','Mediandnds','Meandnds','Genome_Size_Mbp']]
paper_df = paper_df[paper_df.Genus != 'Nonomuraea'].set_index('Genus')

# Create subplots
scale_factor=2
fig, axes = plt.subplots(2, 2,sharey=True,sharex=True, figsize=(8, 6))
fs=10
fs_labels=fs-3

for i in range(2):
    for j in range(2):
        axes[i,j].set_xlim(-0.2,1.1)
        axes[i,j].axvline(x=1, color='black', linestyle='--', linewidth=1)
        axes[i,j].grid(False)
        axes[i,j].set_xticks([0, 1])
        axes[i,j].set_xlabel('FMH dN/dS estimates',fontsize=fs)

axes[0,0].set_title("t=0.02,k=7",fontsize=fs)
axes[0,1].set_title("t=0.05,k=7",fontsize=fs)
axes[1,0].set_title("t=0.1,k=7",fontsize=fs)
axes[1,1].set_title("t=0.2,k=7",fontsize=fs)

# Create scatter plot for containment threshold 0.02
results = pd.concat([data_df1,paper_df],axis=1).dropna().reset_index()
axes[0,0].scatter(results['median1'], results['Genome_Size_Mbp'], color='lightgrey', 
                edgecolor='lightgrey')  # Adjust bins as needed
axes[0,0].set_ylabel('Genome Size (Mbp)',fontsize=fs)

labels=['Buchnera','Blattabacterium','Myxococcus','Actinomyces','Prochlorococcus','Pelagibacter']

for genus in labels:
    coordinates_col3 = results.loc[results['Genus'] == genus, 'median1']
    if not coordinates_col3.empty:
        value_col3 = coordinates_col3.values[0]
    if coordinates_col3.empty:
        value_col3 = None

    coordinates_col1 = results.loc[results['Genus'] == genus, 'Genome_Size_Mbp']
    if not coordinates_col1.empty:
        value_col1 = coordinates_col1.values[0] 
    if coordinates_col1.empty:
        value_col1 = None
    
    #highlight genus of interest
    axes[0,0].scatter(value_col3, value_col1, color='blue',label=genus)
    if value_col1 != None and value_col3 != None:
        if genus == "Blattabacterium" or genus == "Pelagibacter":
            axes[0,0].annotate(genus, (value_col3, value_col1), textcoords="offset points", xytext=(5,-3), ha='left',fontsize=fs_labels)
        elif genus == "Buchnera":
            axes[0,0].annotate(genus, (value_col3, value_col1), textcoords="offset points", xytext=(-5,-3), ha='right',fontsize=fs_labels)
        else:
            axes[0,0].annotate(genus, (value_col3, value_col1), textcoords="offset points", xytext=(5,-3), ha='left',fontsize=fs_labels)

# Create scatter plot for containment threshold 0.05
results = pd.concat([data_df2,paper_df],axis=1).dropna().reset_index()
axes[0,1].scatter(results['median2'], results['Genome_Size_Mbp'], color='lightgrey', 
                edgecolor='lightgrey')  # Adjust bins as needed

for genus in labels:
    coordinates_col3 = results.loc[results['Genus'] == genus, f'median2']
    if not coordinates_col3.empty:
        value_col3 = coordinates_col3.values[0]
    if coordinates_col3.empty:
        value_col3 = None

    coordinates_col1 = results.loc[results['Genus'] == genus, 'Genome_Size_Mbp']
    if not coordinates_col1.empty:
        value_col1 = coordinates_col1.values[0] 
    if coordinates_col1.empty:
        value_col1 = None

    #highlight genus of interest
    axes[0,1].scatter(value_col3, value_col1, color='blue')
    if value_col1 != None and value_col3 != None:
        if genus == "Blattabacterium":
            axes[0,1].annotate(genus, (value_col3, value_col1), textcoords="offset points", xytext=(5,-3), ha='left',fontsize=fs_labels)
        elif genus == "Buchnera":
            axes[0,1].annotate(genus, (value_col3, value_col1), textcoords="offset points", xytext=(-5,-3), ha='right',fontsize=fs_labels)
        else:
            axes[0,1].annotate(genus, (value_col3, value_col1), textcoords="offset points", xytext=(5,-3), ha='left',fontsize=fs_labels) 

# Create scatter plot for containment threshold 0.1
results = pd.concat([data_df3,paper_df],axis=1).dropna().reset_index()
axes[1,0].scatter(results['median3'], results['Genome_Size_Mbp'], color='lightgrey', 
                edgecolor='lightgrey')  # Adjust bins as needed
axes[1,0].set_ylabel('Genome Size (Mbp)',fontsize=fs)

for genus in labels:
    coordinates_col3 = results.loc[results['Genus'] == genus, f'median3']
    print("################################")
    print("Median dN/dS", genus, coordinates_col3)
    if not coordinates_col3.empty:
        value_col3 = coordinates_col3.values[0]
    if coordinates_col3.empty:
        value_col3=None

    coordinates_col1 = results.loc[results['Genus'] == genus, 'Genome_Size_Mbp']
    print("Genome size", genus, coordinates_col1)
    if not coordinates_col1.empty:
        value_col1 = coordinates_col1.values[0]
    if coordinates_col1.empty:
        value_col1=None
    
    #highlight genus of interest
    axes[1,0].scatter(value_col3, value_col1, color='blue')
    if value_col1 != None and value_col3 != None:
        if genus == "Blattabacterium":
            axes[1,0].annotate(genus, (value_col3, value_col1), textcoords="offset points", xytext=(5,-3), ha='left',fontsize=fs_labels)
        else:
            axes[1,0].annotate(genus, (value_col3, value_col1), textcoords="offset points", xytext=(5,-3), ha='left',fontsize=fs_labels) 

# Create scatter plot for containment threshold 0.2
results = pd.concat([data_df4,paper_df],axis=1).dropna().reset_index()
axes[1,1].scatter(results['median3'], results['Genome_Size_Mbp'], color='lightgrey', 
                edgecolor='lightgrey')  # Adjust bins as needed

for genus in labels:
    coordinates_col3 = results.loc[results['Genus'] == genus, f'median3']
    print("################################")
    print("Median dN/dS", genus, coordinates_col3)
    if not coordinates_col3.empty:
        value_col3 = coordinates_col3.values[0]
    if coordinates_col3.empty:
        value_col3=None

    coordinates_col1 = results.loc[results['Genus'] == genus, 'Genome_Size_Mbp']
    print("Genome size", genus, coordinates_col1)
    if not coordinates_col1.empty:
        value_col1 = coordinates_col1.values[0]
    if coordinates_col1.empty:
        value_col1=None
    
    #highlight genus of interest
    axes[1,1].scatter(value_col3, value_col1, color='blue')
    if value_col1 != None and value_col3 != None:
        if genus == "Blattabacterium":
            axes[1,1].annotate(genus, (value_col3, value_col1), textcoords="offset points", xytext=(5,-3), ha='left',fontsize=fs_labels)
        else:
            axes[1,1].annotate(genus, (value_col3, value_col1), textcoords="offset points", xytext=(5,-3), ha='left',fontsize=fs_labels) 


plt.subplots_adjust(wspace=0.01)
plt.tight_layout()
fig.figure.savefig(f"/data/jzr5814/sourmash_dnds_estimation/thesis_figures/figure3b.pdf",bbox_inches='tight') 
