import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
from matplotlib import rcParams, gridspec

# =============================================================================
# PUBLICATION-QUALITY SETTINGS
# =============================================================================

rcParams.update({
    # --- Fonts ---
    'font.family':        'serif',
    'font.size':           8,
    'axes.titlesize':      8,
    'axes.labelsize':      8,
    'xtick.labelsize':     7,
    'ytick.labelsize':     7,
    # --- Lines / axes ---
    'axes.linewidth':      0.6,
    'xtick.major.width':   0.6,
    'ytick.major.width':   0.6,
    'xtick.major.size':    3.0,
    'ytick.major.size':    3.0,
    'lines.linewidth':     0.8,
    # --- Spines ---
    'axes.spines.top':     False,
    'axes.spines.right':   False,
    # --- PDF font embedding ---
    'pdf.fonttype':        42,
    'ps.fonttype':         42,
    'savefig.dpi':         300,
})

# ─── Colors ───────────────────────────────────────────────────────────────────
COL_MAIN = 'blue'
COL_REF  = '#AAAAAA'    # reference lines / gridlines

# ─── Shared labels ────────────────────────────────────────────────────────────
fmh_label   = r"$\mathrm{FMH}\ d_{\mathrm{N}}/d_{\mathrm{S}}$"
codeml_label = r"$\mathrm{CodeML}\ d_{\mathrm{N}}/d_{\mathrm{S}}$"

# =============================================================================
# PARAMETERS
# =============================================================================

ksize      = 7
scaled     = 500
thresholds = [0.02, 0.05, 0.1, 0.2]

PdS_threshold = 0.05
PdN_threshold = 0.23

paper_path = (
    '/data/jzr5814/sourmash_dnds_estimation/tests/data/GTDB/'
    'martinez-gutierrez_2022paper/pgen.1010220.s002.tsv'
)

# =============================================================================
# DATA LOADERS
# =============================================================================

def load_paper_data():
    df = pd.read_csv(paper_path, sep='\t')
    genome_df = df[['Genome', 'Genus', 'Mediandnds', 'Meandnds']].set_index('Genome')
    genus_df  = df[['Genus', 'Mediandnds', 'Meandnds', 'Genome_Size_Mbp']]
    genus_df  = genus_df[genus_df.Genus != 'Nonomuraea'].set_index('Genus')
    return genome_df, genus_df


def load_panelA_data():
    data = {}
    for t in thresholds:
        filepath = (
            f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/'
            f'gtdb_protein_rep_pairwise_{scaled}scale_{t}threshold_k{ksize}/'
            f'fmh_omega_{ksize}_by_genome_group_dSminthresh{PdS_threshold}'
            f'_and_dNminthresh{PdN_threshold}.csv'
        )
        df = pd.read_csv(filepath)[['mean', 'median', 'genome']]
        df['genome'] = df['genome'].apply(
            lambda x: '_'.join(x.split('_')[1:]) if '_' in x else x
        )
        df = df.set_index('genome')
        data[t] = df
    return data


def load_panelB_data():
    data = {}
    for i, t in enumerate(thresholds):
        filepath = (
            f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/'
            f'gtdb_protein_rep_pairwise_{scaled}scale_{t}threshold_k{ksize}/'
            f'fmh_omega_{ksize}_by_genus_group_dSminthresh{PdS_threshold}'
            f'_and_dNminthresh{PdN_threshold}.csv'
        )
        df = pd.read_csv(filepath)[['median', 'genus']]
        df = df.rename(columns={'median': f'median{i}', 'genus': 'Genus'})
        df = df.set_index('Genus')
        data[t] = df
    return data

# =============================================================================
# PANEL A
# =============================================================================

def plot_panelA(axes, data, paper_df):
    xticks = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

    for i, t in enumerate(thresholds):
        ax      = axes[i]
        results = pd.concat([data[t], paper_df], axis=1).dropna()
        corr, _ = scipy.stats.pearsonr(results['median'], results['Mediandnds'])

        ax.scatter(
            results['median'],
            results['Mediandnds'],
            s=8,
            color=COL_MAIN,
            edgecolors='none',
            alpha=0.5,
            rasterized=True,
        )

        ax.set_title(fr"$t={t},\ k_{{\mathrm{{aa}}}}={ksize}$", pad=4)
        ax.set_xlim(-0.04, 0.54)
        ax.set_ylim(-0.01, 0.51)
        ax.set_xticks(xticks)
        ax.set_xticklabels([f'{x:.1f}' for x in xticks])
        ax.tick_params(axis='both', which='both', direction='out')

        # Grid on both axes (vertical line at 0.1 intervals naturally visible)
        ax.yaxis.grid(True, lw=0.25, color='#DDDDDD', zorder=0)
        ax.xaxis.grid(True, lw=0.25, color='#DDDDDD', zorder=0)
        ax.set_axisbelow(True)

        # Pearson r annotation in axis-relative coordinates
        ax.text(
            0.97, 0.04,
            f'$r = {corr:.3f}$',
            transform=ax.transAxes,
            ha='right', va='bottom', fontsize=6.5,
        )

        ax.set_xlabel(fmh_label, labelpad=2)
        ax.tick_params(axis='y', labelleft=(i == 0))

    axes[0].set_ylabel(codeml_label, labelpad=2)

# =============================================================================
# PANEL B
# =============================================================================

def plot_panelB(axes, data, paper_df):
    highlight = [
        'Buchnera', 'Blattabacterium', 'Myxococcus',
        'Actinomyces', 'Prochlorococcus', 'Pelagibacter',
    ]
    label_offsets = {
        'Blattabacterium': ((5, -3),  'left'),
        'Pelagibacter':    ((5, -3),  'left'),
        'Buchnera':        ((-5, -3), 'right'),
    }

    for i, t in enumerate(thresholds):
        ax         = axes[i]
        results    = pd.concat([data[t], paper_df], axis=1).dropna().reset_index()
        median_col = f'median{i}'

        ax.scatter(
            results[median_col],
            results['Genome_Size_Mbp'],
            s=6,
            color=COL_MAIN,
            edgecolors='none',
            alpha=0.05,
            linewidths=0,
            rasterized=True,
        )

        ax.set_xlim(-0.08, 1.08)
        ax.set_title(fr"$t={t},\ k_{{\mathrm{{aa}}}}={ksize}$", pad=4)
        ax.set_xlabel(fmh_label, labelpad=2)
        ax.tick_params(axis='both', which='both', direction='out')

        # Reference line at dN/dS = 1 via grid (1 falls on a natural tick)
        ax.xaxis.grid(True, lw=0.25, color='#DDDDDD', zorder=0)
        ax.yaxis.grid(True, lw=0.25, color='#DDDDDD', zorder=0)
        ax.set_axisbelow(True)
        # Explicit reference line at dN/dS = 1 (thicker so it reads clearly)
        ax.axvline(x=1, lw=0.8, ls='--', color=COL_REF, zorder=2)

        ax.tick_params(axis='y', labelleft=(i == 0))

        # Highlighted genera
        for genus in highlight:
            row = results[results['Genus'] == genus]
            if row.empty:
                continue
            x = row[median_col].values[0]
            y = row['Genome_Size_Mbp'].values[0]
            ax.scatter(x, y, s=20, color=COL_MAIN, edgecolors='none', zorder=3)
            (dx, dy), ha = label_offsets.get(genus, ((5, -3), 'left'))
            ax.annotate(
                genus, (x, y),
                xytext=(dx, dy),
                textcoords='offset points',
                ha=ha, fontsize=6,
            )

    axes[0].set_ylabel('Genome Size (Mbp)', labelpad=2)

# =============================================================================
# MAIN FIGURE
# =============================================================================

def make_figure():
    genome_df, genus_df = load_paper_data()
    panelA_data = load_panelA_data()
    panelB_data = load_panelB_data()

    fig    = plt.figure(figsize=(7.5, 5.5))
    outer  = gridspec.GridSpec(2, 1, hspace=0.50)

    gsA   = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=outer[0], wspace=0.10)
    axesA = [plt.subplot(gsA[i]) for i in range(4)]
    plot_panelA(axesA, panelA_data, genome_df)

    gsB   = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=outer[1], wspace=0.10)
    axesB = [plt.subplot(gsB[i]) for i in range(4)]
    plot_panelB(axesB, panelB_data, genus_df)

    plt.subplots_adjust(left=0.09, right=0.98, top=0.95, bottom=0.10)

    fig.text(0.01, 0.97, 'A', fontsize=10, fontweight='bold', va='top')
    fig.text(0.01, 0.48, 'B', fontsize=10, fontweight='bold', va='top')

    fig.savefig(
        '/data/jzr5814/repositories/dnds_using_fmh_reproducibles/'
        'manuscript_figures/updated_pdf/figure4_now5.pdf',
        bbox_inches='tight',
        dpi=300,
    )


if __name__ == '__main__':
    make_figure()