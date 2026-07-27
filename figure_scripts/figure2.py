import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
from matplotlib import rcParams
import matplotlib.pyplot as plt

rcParams['font.family'] = 'serif'
#ylabeling = r"$\mathrm{FMH}\ d_{\mathrm{N}}/d_{\mathrm{S}}$"
ylabeling = "FMH Omega"
# =========================
# -------- HELPERS --------
# =========================

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
    pos = pd.concat([fmh_pos, ng_pos], axis=1).dropna()
    neg = pd.concat([fmh_neg, ng_neg], axis=1).dropna()
    combined = pd.concat([pos, neg], ignore_index=True)
    r, _ = scipy.stats.pearsonr(combined['Ka/Ks'], combined['dN/dS'])
    return pos, neg, r

def plot_panel(ax, pos, neg, r, fs=12):
    ax.scatter(neg['dN/dS'], neg['Ka/Ks'], color='blue', alpha=0.5)
    ax.scatter(pos['dN/dS'], pos['Ka/Ks'], color='darkorange', alpha=0.5)
    tmp = [neg['Ka/Ks'].min(), pos['Ka/Ks'].max()]
    ax.plot(tmp, tmp, linestyle='--', color='grey')
    ax.text(4, 9, f'Pearson R: {r:.3f}', fontsize=fs-2)

def format_ax(ax, xlim, ylim, fs=12, xlabel=False, ylabel=False):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xticks([0,1,2,5,10])
    ax.set_yticks([0,1,2,5,10])
    ax.tick_params(axis='both', labelsize=fs-2)
    if xlabel:
        ax.set_xlabel(ylabeling, fontsize=fs)
    if ylabel:
        ax.set_ylabel('NG86', fontsize=fs)

# =========================
# -------- FILES A --------
# =========================

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

# =========================
# -------- FILES B --------
# =========================

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

# =========================
# -------- FIGURE ---------
# =========================

# 1. Use a square figure size
fig = plt.figure(figsize=(12, 12))

# 2. Set height_ratios so all data rows (0,1,2 and 4) are EXACTLY 1.
# The spacer (row 3) is 0.4.
# hspace=0.1 provides a small uniform gap between all rows.
gs = fig.add_gridspec(
    7,
    4,
    height_ratios=[1,1,1,0.3,1,0.3,1],
    hspace=0.4,
    wspace=0.1
)

# --- Define axes_A (Rows 0, 1, 2) ---
ax_master_A = fig.add_subplot(gs[0, 0])
axes_A = [[ax_master_A if (i == 0 and j == 0) 
           else fig.add_subplot(gs[i, j], sharex=ax_master_A, sharey=ax_master_A) 
           for j in range(4)] for i in range(3)]

# --- Define axes_B (Row 4, Columns 0, 1, 2) ---
# We MUST use the 4-column grid here to keep widths consistent with Panel A
ax_master_B = fig.add_subplot(gs[4, 0])
axes_B = [ax_master_B if j == 0 
          else fig.add_subplot(gs[4, j], sharex=ax_master_B, sharey=ax_master_B) 
          for j in range(3)]

# Panel C axes
ax_master_C = fig.add_subplot(gs[6,0])
axes_C = [
    ax_master_C if i == 0
    else fig.add_subplot(
        gs[6,i],
        sharex=ax_master_C,
        sharey=ax_master_C
    )
    for i in range(3)
]

ghost_ax_C = fig.add_subplot(gs[5,3])
ghost_ax_C.set_visible(False)

# --- Add the Ghost Plot in Column 3 ---
# This "holds" the 4th column open so plots 0, 1, and 2 don't stretch to fill the row
ghost_ax = fig.add_subplot(gs[4, 3])
ghost_ax.set_visible(False)

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
     r"$k_{\mathrm{aa}} = 21,\ 20001\ \mathrm{bp}$"]
]

titles_B = [
    r"$k_{\mathrm{aa}} = 7$",
    r"$k_{\mathrm{aa}} = 9$",
    r"$k_{\mathrm{aa}} = 11$"
]


# =========================
# -------- PANEL A --------
# =========================

for i in range(3):        # i = row (0, 1, 2)
    for j in range(4):    # j = col (0, 1, 2, 3)
        ax = axes_A[i][j]

        ax.set_title(titles_A[i][j], fontsize=14, pad=8)
        
        # Load data (using your existing logic)
        fmh_pos = load_fmh(fmh_dnds_positive_files[j][i], 'ref_gene')
        fmh_neg = load_fmh(fmh_dnds_negative_files[j][i], 'ref_gene')
        ng_pos = load_ng(ng_dnds_positive_files[i])
        ng_neg = load_ng(ng_dnds_negative_files[i])

        pos, neg, r = combine_and_corr(fmh_pos, fmh_neg, ng_pos, ng_neg)

        # Plot the data
        plot_panel(ax, pos, neg, r)
        
        # Formatting Logic:
        # Show ylabel only if leftmost column (j == 0)
        # Show xlabel only if bottom row (i == 2)
        is_left_edge = (j == 0)
        is_bottom_edge = (i == 2)
        
        format_ax(ax, (0, 11), (0, 11),
                  xlabel=is_bottom_edge,
                  ylabel=is_left_edge)

        # Explicitly hide tick labels for inner plots to keep it clean
        if not is_left_edge:
            ax.tick_params(labelleft=False)
        if not is_bottom_edge:
            ax.tick_params(labelbottom=False)

# =========================
# -------- PANEL B --------
# =========================

for k in range(3):
    ax = axes_B[k]

    ax.set_title(titles_B[k], fontsize=14, pad=10)

    # Load data (using your existing logic)
    fmh_pos = load_fmh(fmh_dnds_positive_files_B[k], 'positive_0.01_ref', upper=8)
    fmh_neg = load_fmh(fmh_dnds_negative_files_B[k], 'negative_0.01_ref', upper=8)

    ng_pos = load_ng(ng86_dnds_positive_files, upper=8, ref_label='positive_0.01_ref')
    ng_neg = load_ng(ng86_dnds_negative_files, upper=8, ref_label='negative_0.01_ref')

    pos, neg, r = combine_and_corr(fmh_pos, fmh_neg, ng_pos, ng_neg)

    # Plot the data
    plot_panel(ax, pos, neg, r)
    
    # Formatting Logic for Row B:
    # ylabel only for the first one (k == 0)
    # xlabel for all of them (since they are the bottom-most row of the figure)
    is_left_edge = (k == 0)
    
    format_ax(ax, (0, 11), (0, 11),
              xlabel=True, 
              ylabel=is_left_edge)

    # Explicitly hide Y-tick labels for the 2nd and 3rd plots
    if not is_left_edge:
        ax.tick_params(labelleft=False)

# =========================
# -------- PANEL C --------
# =========================

def load_hist_fmh(file, ref_label='ref_ecoli', upper=10):
    df = pd.read_csv(file, sep=',')
    df = df[df['query_name'] == ref_label][['match_name', 'dN/dS']]
    df = df[(df['dN/dS'] >= 0) & (df['dN/dS'] <= upper)]
    return df


fmh_negative_files_C = [
    "/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k9_10scale_cores1_includes_ref/fmh_omega_9.csv",
    "/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k9_100scale_cores1_includes_ref/fmh_omega_9.csv",
    "/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k9_1000scale_cores1_includes_ref/fmh_omega_9.csv"
]

fmh_positive_files_C = [
    "/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k9_10scale_cores1_includes_ref/fmh_omega_9.csv",
    "/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k9_100scale_cores1_includes_ref/fmh_omega_9.csv",
    "/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k9_1000scale_cores1_includes_ref/fmh_omega_9.csv"
]


negative_C = [
    load_hist_fmh(f)
    for f in fmh_negative_files_C
]

positive_C = [
    load_hist_fmh(f)
    for f in fmh_positive_files_C
]

scales=[10,100,1000]

for i in range(3):

    ax = axes_C[i]

    ax.hist(
        negative_C[i]['dN/dS'],
        bins=50,
        range=(0,8),
        alpha=0.5,
        color='blue',
        edgecolor='blue'
    )

    ax.hist(
        positive_C[i]['dN/dS'],
        bins=50,
        range=(0,8),
        alpha=0.5,
        color='darkorange',
        edgecolor='darkorange'
    )

    ax.axvline(
        x=1,
        linestyle='--',
        color='grey'
    )

    ax.set_title(
        f"scaled={scales[i]}",
        fontsize=14
    )

    ax.set_xlim(-0.5,8.5)
    ax.set_xticks([0,2,4,6,8])
    ax.tick_params(axis='both', labelsize=10)

    ax.set_xlabel(
        ylabeling,
        fontsize=12
    )

    if i != 0:
        ax.tick_params(axis='y', labelleft=False)

    if i == 0:
        ax.set_ylabel(
            "Frequency",
            fontsize=12
        )
    

# =========================
# -------- FINAL ---------
# =========================

plt.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.08)

fig.text(0.02, 0.98, 'A', fontsize=20, fontweight='bold')
fig.text(0.02, 0.46, 'B', fontsize=20, fontweight='bold')
fig.text(0.02, 0.22, "C", fontsize=20, fontweight="bold")

fig.savefig("/data/jzr5814/repositories/dnds_using_fmh_reproducibles/manuscript_figures/updated_pdf/figure2.png",
            bbox_inches='tight')
