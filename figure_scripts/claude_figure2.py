import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
from matplotlib import rcParams
import matplotlib.patches as mpatches

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
    # --- Spines (remove top/right globally) ---
    'axes.spines.top':     False,
    'axes.spines.right':   False,
    # --- PDF font embedding (required by most journals) ---
    'pdf.fonttype':        42,
    'ps.fonttype':         42,
    'savefig.dpi':         300,
})

# ─── Color palette (blue / amber – colorblind-friendly) ───────────────────────
COL_NEG  = 'blue'    # negative selection
COL_POS  = 'darkorange'    # positive selection
COL_REF  = '#AAAAAA'    # 1:1 diagonal / dN/dS = 1 reference

# ─── Shared axis labels ───────────────────────────────────────────────────────
fmh_label = r"$\mathrm{FMH}\ d_{\mathrm{N}}/d_{\mathrm{S}}$"
ng_label  = r"KaKs Calculator (NG86)"

# =============================================================================
# HELPERS
# =============================================================================

def load_fmh(file, ref_label, upper=10):
    df = pd.read_csv(file, sep=',')
    df = df[df['A'] == ref_label][['B', 'dN/dS']]
    df = df[(df['dN/dS'] > 0) & (df['dN/dS'] < upper)]
    return df

def load_ng(file, upper=10, ref_label=None):
    df = pd.read_csv(file, sep='\t')[['Sequence', 'Ka/Ks']]
    if ref_label:
        df = df[df['Sequence'] == ref_label]
    df = df[(df['Ka/Ks'] > 0) & (df['Ka/Ks'] < upper)]
    return df

def combine_and_corr(fmh_pos, fmh_neg, ng_pos, ng_neg):
    pos     = pd.concat([fmh_pos, ng_pos], axis=1).dropna()
    neg     = pd.concat([fmh_neg, ng_neg], axis=1).dropna()
    combined = pd.concat([pos, neg], ignore_index=True)
    r, _    = scipy.stats.pearsonr(combined['Ka/Ks'], combined['dN/dS'])
    return pos, neg, r

def plot_scatter(ax, pos, neg, r):
    """Draw one scatter panel: negative (blue) and positive (amber)."""
    ax.scatter(
        neg['dN/dS'], neg['Ka/Ks'],
        s=5, color=COL_NEG, alpha=0.55,
        edgecolors='none', linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        pos['dN/dS'], pos['Ka/Ks'],
        s=5, color=COL_POS, alpha=0.55,
        edgecolors='none', linewidths=0,
        rasterized=True,
    )
    ax.plot([0, 11], [0, 11], ls='--', lw=0.7, color=COL_REF, zorder=0)
    # Annotation in axis-relative coords → robust to any xlim/ylim
    ax.text(
        0.97, 0.04,
        f'$r = {r:.3f}$',
        transform=ax.transAxes,
        ha='right', va='bottom', fontsize=6.5,
    )

def format_scatter(ax, xlim=(0, 11), ylim=(0, 11),
                   xlabel=False, ylabel=False):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xticks([0, 1, 2, 5, 10])
    ax.set_yticks([0, 1, 2, 5, 10])
    ax.tick_params(axis='both', which='both', direction='out')
    # Subtle horizontal grid only
    ax.yaxis.grid(True, lw=0.25, color='#DDDDDD', zorder=0)
    ax.xaxis.grid(True, lw=0.25, color='#DDDDDD', zorder=0)
    ax.set_axisbelow(True)
    if xlabel:
        ax.set_xlabel(fmh_label, labelpad=2)
    if ylabel:
        ax.set_ylabel(ng_label, labelpad=2)

# =============================================================================
# FILES – PANEL A
# =============================================================================

fmh_dnds_k5_negative_files = [
    "../data/figure2_datasets/fmh_omega_5_negative_5001bp.csv",
    "../data/figure2_datasets/fmh_omega_5_negative_10002bp.csv",
    "../data/figure2_datasets/fmh_omega_5_negative_20001bp.csv",
]
fmh_dnds_k7_negative_files = [
    "../data/figure2_datasets/fmh_omega_7_negative_5001bp.csv",
    "../data/figure2_datasets/fmh_omega_7_negative_10002bp.csv",
    "../data/figure2_datasets/fmh_omega_7_negative_20001bp.csv",
]
fmh_dnds_k15_negative_files = [
    "../data/figure2_datasets/fmh_omega_15_negative_5001bp.csv",
    "../data/figure2_datasets/fmh_omega_15_negative_10002bp.csv",
    "../data/figure2_datasets/fmh_omega_15_negative_20001bp.csv",
]
fmh_dnds_k21_negative_files = [
    "../data/figure2_datasets/fmh_omega_21_negative_5001bp.csv",
    "../data/figure2_datasets/fmh_omega_21_negative_10002bp.csv",
    "../data/figure2_datasets/fmh_omega_21_negative_20001bp.csv",
]
fmh_dnds_k5_positive_files = [
    "../data/figure2_datasets/fmh_omega_5_positive_5001bp.csv",
    "../data/figure2_datasets/fmh_omega_5_positive_10002bp.csv",
    "../data/figure2_datasets/fmh_omega_5_positive_20001bp.csv",
]
fmh_dnds_k7_positive_files = [
    "../data/figure2_datasets/fmh_omega_7_positive_5001bp.csv",
    "../data/figure2_datasets/fmh_omega_7_positive_10002bp.csv",
    "../data/figure2_datasets/fmh_omega_7_positive_20001bp.csv",
]
fmh_dnds_k15_positive_files = [
    "../data/figure2_datasets/fmh_omega_15_positive_5001bp.csv",
    "../data/figure2_datasets/fmh_omega_15_positive_10002bp.csv",
    "../data/figure2_datasets/fmh_omega_15_positive_20001bp.csv",
]
fmh_dnds_k21_positive_files = [
    "../data/figure2_datasets/fmh_omega_21_positive_5001bp.csv",
    "../data/figure2_datasets/fmh_omega_21_positive_10002bp.csv",
    "../data/figure2_datasets/fmh_omega_21_positive_20001bp.csv",
]

fmh_dnds_negative_files = [
    fmh_dnds_k5_negative_files,
    fmh_dnds_k7_negative_files,
    fmh_dnds_k15_negative_files,
    fmh_dnds_k21_negative_files,
]
fmh_dnds_positive_files = [
    fmh_dnds_k5_positive_files,
    fmh_dnds_k7_positive_files,
    fmh_dnds_k15_positive_files,
    fmh_dnds_k21_positive_files,
]

ng_dnds_positive_files = [
    "../data/figure2_datasets/positive_selection_queries_5001_0.01.axt.kaks",
    "../data/figure2_datasets/positive_selection_queries_10002_0.01.axt.kaks",
    "../data/figure2_datasets/positive_selection_queries_20001_0.01.axt.kaks",
]
ng_dnds_negative_files = [
    "../data/figure2_datasets/negative_selection_queries_5001_0.01.axt.kaks",
    "../data/figure2_datasets/negative_selection_queries_10002_0.01.axt.kaks",
    "../data/figure2_datasets/negative_selection_queries_20001_0.01.axt.kaks",
]

# =============================================================================
# FILES – PANEL B
# =============================================================================

fmh_dnds_positive_files_B = [
    "../data/figure2_datasets/fmh_omega_7_positive_LAMA3.csv",
    "../data/figure2_datasets/fmh_omega_9_positive_LAMA3.csv",
    "../data/figure2_datasets/fmh_omega_11_positive_LAMA3.csv",
]
fmh_dnds_negative_files_B = [
    "../data/figure2_datasets/fmh_omega_7_negative_LAMA3.csv",
    "../data/figure2_datasets/fmh_omega_9_negative_LAMA3.csv",
    "../data/figure2_datasets/fmh_omega_11_negative_LAMA3.csv",
]
ng86_dnds_positive_files = "../data/figure2_datasets/kaks_sequences_positive_LAMA3.axt.kaks"
ng86_dnds_negative_files = "../data/figure2_datasets/kaks_sequences_negative_LAMA3.axt.kaks"

# =============================================================================
# FILES – PANEL C
# =============================================================================

def load_hist_fmh(file, ref_label='ref_ecoli', upper=10):
    df = pd.read_csv(file, sep=',')
    df = df[df['query_name'] == ref_label][['match_name', 'dN/dS']]
    df = df[(df['dN/dS'] >= 0) & (df['dN/dS'] <= upper)]
    return df

fmh_negative_files_C = [
    "/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k9_10scale_cores1_includes_ref/fmh_omega_9.csv",
    "/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k9_100scale_cores1_includes_ref/fmh_omega_9.csv",
    "/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k9_1000scale_cores1_includes_ref/fmh_omega_9.csv",
]
fmh_positive_files_C = [
    "/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k9_10scale_cores1_includes_ref/fmh_omega_9.csv",
    "/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k9_100scale_cores1_includes_ref/fmh_omega_9.csv",
    "/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k9_1000scale_cores1_includes_ref/fmh_omega_9.csv",
]

negative_C = [load_hist_fmh(f) for f in fmh_negative_files_C]
positive_C = [load_hist_fmh(f) for f in fmh_positive_files_C]

# =============================================================================
# FIGURE LAYOUT
# =============================================================================

PANEL_FS = 10   # bold A / B / C label size

fig = plt.figure(figsize=(7.5, 9.5))

# 7 rows: A×3 | spacer | B | spacer | C
gs = fig.add_gridspec(
    7, 4,
    height_ratios=[1, 1, 1, 0.22, 1, 0.22, 1],
    hspace=0.20,
    wspace=0.08,
)

# ── Panel A (rows 0-2, all 4 cols) ──────────────────────────────────────────
ax_master_A = fig.add_subplot(gs[0, 0])
axes_A = [
    [ax_master_A if (i == 0 and j == 0)
     else fig.add_subplot(gs[i, j], sharex=ax_master_A, sharey=ax_master_A)
     for j in range(4)]
    for i in range(3)
]

# ── Panel B (row 4, cols 0-2; col 3 = legend) ───────────────────────────────
ax_master_B = fig.add_subplot(gs[4, 0])
axes_B = [
    ax_master_B if j == 0
    else fig.add_subplot(gs[4, j], sharex=ax_master_B, sharey=ax_master_B)
    for j in range(3)
]
legend_ax = fig.add_subplot(gs[4, 3])
legend_ax.axis('off')

# ── Panel C (row 6, cols 0-2; col 3 = ghost) ────────────────────────────────
ax_master_C = fig.add_subplot(gs[6, 0])
axes_C = [
    ax_master_C if i == 0
    else fig.add_subplot(gs[6, i], sharex=ax_master_C, sharey=ax_master_C)
    for i in range(3)
]
ghost_C = fig.add_subplot(gs[6, 3])
ghost_C.axis('off')

# ── Titles ───────────────────────────────────────────────────────────────────
titles_A = [
    [r"$k_{\mathrm{aa}} = 5,\ 5001\ \mathrm{bp}$",
     r"$k_{\mathrm{aa}} = 7,\ 5001\ \mathrm{bp}$",
     r"$k_{\mathrm{aa}} = 15,\ 5001\ \mathrm{bp}$",
     r"$k_{\mathrm{aa}} = 21,\ 5001\ \mathrm{bp}$"],
    [r"$k_{\mathrm{aa}} = 5,\ 10002\ \mathrm{bp}$",
     r"$k_{\mathrm{aa}} = 7,\ 10002\ \mathrm{bp}$",
     r"$k_{\mathrm{aa}} = 15,\ 10002\ \mathrm{bp}$",
     r"$k_{\mathrm{aa}} = 21,\ 10002\ \mathrm{bp}$"],
    [r"$k_{\mathrm{aa}} = 5,\ 20001\ \mathrm{bp}$",
     r"$k_{\mathrm{aa}} = 7,\ 20001\ \mathrm{bp}$",
     r"$k_{\mathrm{aa}} = 15,\ 20001\ \mathrm{bp}$",
     r"$k_{\mathrm{aa}} = 21,\ 20001\ \mathrm{bp}$"],
]
titles_B = [
    r"$k_{\mathrm{aa}} = 7$",
    r"$k_{\mathrm{aa}} = 9$",
    r"$k_{\mathrm{aa}} = 11$",
]
scales_C = [10, 100, 1000]

# =============================================================================
# PANEL A
# =============================================================================

for i in range(3):
    for j in range(4):
        ax = axes_A[i][j]
        ax.set_title(titles_A[i][j], pad=4)

        fmh_pos = load_fmh(fmh_dnds_positive_files[j][i], 'ref_gene')
        fmh_neg = load_fmh(fmh_dnds_negative_files[j][i], 'ref_gene')
        ng_pos  = load_ng(ng_dnds_positive_files[i])
        ng_neg  = load_ng(ng_dnds_negative_files[i])
        pos, neg, r = combine_and_corr(fmh_pos, fmh_neg, ng_pos, ng_neg)

        plot_scatter(ax, pos, neg, r)
        format_scatter(
            ax,
            xlabel=(i == 2),
            ylabel=(j == 0),
        )
        if j != 0:
            ax.tick_params(labelleft=False)
        if i != 2:
            ax.tick_params(labelbottom=False)

# =============================================================================
# PANEL B
# =============================================================================

for k in range(3):
    ax = axes_B[k]
    ax.set_title(titles_B[k], pad=4)

    fmh_pos = load_fmh(fmh_dnds_positive_files_B[k], 'positive_0.01_ref', upper=8)
    fmh_neg = load_fmh(fmh_dnds_negative_files_B[k], 'negative_0.01_ref', upper=8)
    ng_pos  = load_ng(ng86_dnds_positive_files, upper=8, ref_label='positive_0.01_ref')
    ng_neg  = load_ng(ng86_dnds_negative_files, upper=8, ref_label='negative_0.01_ref')
    pos, neg, r = combine_and_corr(fmh_pos, fmh_neg, ng_pos, ng_neg)

    plot_scatter(ax, pos, neg, r)
    format_scatter(ax, xlabel=True, ylabel=(k == 0))
    if k != 0:
        ax.tick_params(labelleft=False)

# ── Legend in the 4th column of Panel B ──────────────────────────────────────
neg_patch = mpatches.Patch(
    facecolor=COL_NEG, alpha=0.65, linewidth=0,
    label='Negative selection',
)
pos_patch = mpatches.Patch(
    facecolor=COL_POS, alpha=0.65, linewidth=0,
    label='Positive selection',
)
legend_ax.legend(
    handles=[neg_patch, pos_patch],
    loc='center left',
    frameon=False,
    fontsize=7.5,
    title='Selection type',
    title_fontsize=7.5,
    handlelength=1.0,
    handleheight=0.9,
    borderpad=0,
)

# =============================================================================
# PANEL C
# =============================================================================

for i in range(3):
    ax = axes_C[i]

    ax.hist(
        negative_C[i]['dN/dS'],
        bins=50, range=(0, 8),
        alpha=0.60, color=COL_NEG, linewidth=0,
    )
    ax.hist(
        positive_C[i]['dN/dS'],
        bins=50, range=(0, 8),
        alpha=0.60, color=COL_POS, linewidth=0,
    )
    ax.axvline(x=1, ls='--', lw=0.8, color=COL_REF, zorder=3)

    ax.set_title(
        r"$\mathrm{scale} = " + str(scales_C[i]) + r"$",
        pad=4,
    )
    ax.set_xlim(-0.5, 8.5)
    ax.set_xticks([0, 2, 4, 6, 8])
    ax.tick_params(axis='both', which='both', direction='out')
    ax.yaxis.grid(True, lw=0.25, color='#DDDDDD', zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlabel(fmh_label, labelpad=2)

    if i == 0:
        ax.set_ylabel('Frequency', labelpad=2)
    else:
        ax.tick_params(labelleft=False)

# =============================================================================
# FINAL TOUCHES
# =============================================================================

plt.subplots_adjust(left=0.10, right=0.98, top=0.97, bottom=0.05)

# Panel labels – positioned left of the plot area, vertically centered
# on each panel group. Using figure-fraction coordinates.
for label, x, y in [('A', 0.01, 0.985), ('B', 0.01, 0.455), ('C', 0.01, 0.215)]:
    fig.text(x, y, label,
             fontsize=PANEL_FS, fontweight='bold',
             va='top', ha='left')

fig.savefig(
    "/data/jzr5814/repositories/dnds_using_fmh_reproducibles/manuscript_figures/updated_pdf/figure2.pdf",
    bbox_inches='tight',
    dpi=300,
)