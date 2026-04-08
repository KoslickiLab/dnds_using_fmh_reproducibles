import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
from matplotlib import rcParams, gridspec

# -----------------------------
# GLOBAL SETTINGS
# -----------------------------
rcParams['font.family'] = 'serif'

ksize = 7
scaled = 500
thresholds = [0.02, 0.05, 0.1, 0.2]

PdS_threshold = 0.05
PdN_threshold = 0.23

paper_path = '/data/jzr5814/sourmash_dnds_estimation/tests/data/GTDB/martinez-gutierrez_2022paper/pgen.1010220.s002.tsv'


# -----------------------------
# LOAD PAPER DATA
# -----------------------------
def load_paper_data():
    df = pd.read_csv(paper_path, sep='\t')
    genome_df = df[['Genome', 'Genus', 'Mediandnds', 'Meandnds']].set_index('Genome')
    genus_df = df[['Genus', 'Mediandnds', 'Meandnds', 'Genome_Size_Mbp']]
    genus_df = genus_df[genus_df.Genus != 'Nonomuraea'].set_index('Genus')
    return genome_df, genus_df


# -----------------------------
# LOAD PANEL A DATA
# -----------------------------
def load_panelA_data():
    data = {}

    for t in thresholds:
        filepath = f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_{scaled}scale_{t}threshold_k{ksize}/fmh_omega_{ksize}_by_genome_group_dSminthresh{PdS_threshold}_and_dNminthresh{PdN_threshold}.csv'

        df = pd.read_csv(filepath)[['mean', 'median', 'genome']]
        df['genome'] = df['genome'].apply(lambda x: '_'.join(x.split('_')[1:]) if '_' in x else x)
        df = df.set_index('genome')

        data[t] = df

    return data


# -----------------------------
# LOAD PANEL B DATA
# -----------------------------
def load_panelB_data():
    data = {}

    for i, t in enumerate(thresholds):
        filepath = f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_{scaled}scale_{t}threshold_k{ksize}/fmh_omega_{ksize}_by_genus_group_dSminthresh{PdS_threshold}_and_dNminthresh{PdN_threshold}.csv'

        df = pd.read_csv(filepath)[['median', 'genus']]
        df = df.rename(columns={'median': f'median{i}', 'genus': 'Genus'})
        df = df.set_index('Genus')

        data[t] = df

    return data


# -----------------------------
# PANEL A PLOTTING
# -----------------------------
def plot_panelA(axes, data, paper_df):
    fs = 14
    xticks = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

    for i, t in enumerate(thresholds):
        results = pd.concat([data[t], paper_df], axis=1).dropna()

        corr, _ = scipy.stats.pearsonr(results['median'], results['Mediandnds'])

        axes[i].scatter(
            results['median'],
            results['Mediandnds'],
            color='blue',
            edgecolor='blue',
            alpha=0.5
        )

        axes[i].set_title(f"t={t}, k={ksize}", fontsize=fs)
        axes[i].set_xlim(-0.04, 0.54)
        axes[i].set_ylim(-0.01, 0.51)

        axes[i].text(0.2, 0.4, f'R = {corr:.3f}', fontsize=fs-2)

        axes[i].tick_params(labelsize=fs-2)

        axes[i].set_xlabel('FMH dN/dS estimates', fontsize=fs)

        axes[i].tick_params(axis='y', labelleft=(i == 0))
        axes[i].set_xticks(xticks)
        axes[i].set_xticklabels([f"{x:.1f}" for x in xticks])

    axes[0].set_ylabel('CodeML dN/dS estimates', fontsize=fs)


# -----------------------------
# PANEL B PLOTTING
# -----------------------------
def plot_panelB(axes, data, paper_df):
    fs = 14
    fs_labels = fs - 4

    highlight = ['Buchnera','Blattabacterium','Myxococcus',
                 'Actinomyces','Prochlorococcus','Pelagibacter']

    for i, t in enumerate(thresholds):
        results = pd.concat([data[t], paper_df], axis=1).dropna().reset_index()

        median_col = f'median{i}'

        axes[i].scatter(
            results[median_col],
            results['Genome_Size_Mbp'],
            color='lightgrey',
            edgecolor='lightgrey'
        )

        axes[i].set_xlim(-0.05, 1.05)
        axes[i].axvline(x=1, linestyle='--', color='black')
        axes[i].set_title(f"t={t}, k={ksize}", fontsize=fs)
        axes[i].set_xlabel('FMH dN/dS estimates', fontsize=fs)
        axes[i].tick_params(axis='y', labelleft=(i == 0))
        

        # highlight genera
        for genus in highlight:
            row = results[results['Genus'] == genus]

            if not row.empty:
                x = row[median_col].values[0]
                y = row['Genome_Size_Mbp'].values[0]

                axes[i].scatter(x, y, color='blue')
                axes[i].annotate(genus, (x, y),
                                 xytext=(5, -3),
                                 textcoords='offset points',
                                 fontsize=fs_labels)

    axes[0].set_ylabel('Genome Size (Mbp)', fontsize=fs)


# -----------------------------
# MAIN FIGURE
# -----------------------------
def make_figure():
    genome_df, genus_df = load_paper_data()
    panelA_data = load_panelA_data()
    panelB_data = load_panelB_data()

    fig = plt.figure(figsize=(16, 8))
    outer = gridspec.GridSpec(2, 1, hspace=0.5)

    # Panel A
    gsA = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=outer[0], wspace=0.1)
    axesA = [plt.subplot(gsA[i]) for i in range(4)]
    plot_panelA(axesA, panelA_data, genome_df)

    # Panel B
    gsB = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=outer[1], wspace=0.1)
    axesB = [plt.subplot(gsB[i]) for i in range(4)]
    plot_panelB(axesB, panelB_data, genus_df)

    # Panel labels
    fig.text(0.07, 0.92, 'A', fontsize=18, fontweight='bold')
    fig.text(0.07, 0.47, 'B', fontsize=18, fontweight='bold')

    plt.savefig("/data/jzr5814/repositories/dnds_using_fmh_reproducibles/manuscript_figures/updated_pdf/figure4_now5.pdf", bbox_inches='tight')
    plt.show()


# -----------------------------
# RUN
# -----------------------------
if __name__ == "__main__":
    make_figure()
