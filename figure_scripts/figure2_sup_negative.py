import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
from matplotlib import rcParams
import matplotlib.pyplot as plt

rcParams['font.family'] = 'serif'

# =========================
# -------- HELPERS --------
# =========================

def load_fmh(file, ref_label, upper=2):
    df = pd.read_csv(file, sep=',')
    df = df[df['A'] == ref_label][['B', 'dN/dS']]
    df = df[(df['dN/dS'] > 0) & (df['dN/dS'] < upper)]
    return df

def load_ng(file, upper=2, ref_label=None):
    df = pd.read_csv(file, sep='\t')[['Sequence', 'Ka/Ks']]
    if ref_label:
        df = df[df['Sequence'] == ref_label]
    df = df[(df['Ka/Ks'] > 0) & (df['Ka/Ks'] < upper)]
    return df

def combine_and_corr(fmh_neg, ng_neg):
    neg = pd.concat([fmh_neg, ng_neg], axis=1).dropna()
    r, _ = scipy.stats.pearsonr(neg['Ka/Ks'], neg['dN/dS'])
    return neg, r

def plot_panel(ax, neg, r, fs=13):
    ax.scatter(neg['dN/dS'], neg['Ka/Ks'], color='blue', alpha=0.5)
    ax.text(0.6, 1.3, f'Pearson R: {r:.3f}', fontsize=fs-2)

def format_ax(ax, xlim, ylim, fs=13, xlabel=False, ylabel=False):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xticks([0,1])
    ax.set_yticks([0,1])
    ax.tick_params(axis='both', labelsize=fs-2)
    if xlabel:
        ax.set_xlabel('FMH dN/dS estimates', fontsize=fs)
    if ylabel:
        ax.set_ylabel('NG86 dN/dS estimates', fontsize=fs)

# =========================
# -------- FILES A --------
# =========================

fmh_dnds_k5_negative_files = [
    '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/negative/k5/fmh_omega_5.csv',
    '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/ksizes/5/negative/fmh_omega_5.csv',
    '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/negative/k5/fmh_omega_5.csv'
]

fmh_dnds_k7_negative_files = [
    '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/negative/fmh_omega_7.csv',
    '/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/fmh_dnds_sketch_protein/negative_selection_redo_sketch_protein_using_faa/dnds.k7.approximations_included.csv',
    '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/negative/fmh_omega_7.csv'
]

fmh_dnds_k15_negative_files = [
    '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/negative/k15/fmh_omega_15.csv',
    '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/ksizes/15/negative/fmh_omega_15.csv',
    '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/negative/k15/fmh_omega_15.csv'
]

fmh_dnds_k21_negative_files = [
    '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/negative/k21/fmh_omega_21.csv',
    '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/ksizes/21/negative/fmh_omega_21.csv',
    '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/negative/k21/fmh_omega_21.csv'
]

fmh_dnds_negative_files = [
    fmh_dnds_k5_negative_files,
    fmh_dnds_k7_negative_files,
    fmh_dnds_k15_negative_files,
    fmh_dnds_k21_negative_files
]

ng_dnds_negative_files = [
    '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/5001/kaks_NG/negative_selection_queries_5001_0.01.axt.kaks',
    '/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/kaks_NG/negative_selection_queries_10002_0.01.axt.kaks',
    '/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/20001/kaks_NG/negative_selection_queries_20001_0.01.axt.kaks'
]

# =========================
# -------- FILES B --------
# =========================


fmh_dnds_negative_files_B = [
    '/data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/negative/fmh_k7/fmh_omega_7_modified.csv',
    '/data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/negative/fmh_k9/fmh_omega_9.csv',
    '/data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/negative/fmh_k11/fmh_omega_11.csv'
]

ng86_dnds_negative_files = '/data/jzr5814/sourmash_dnds_estimation/tests/test/real_sequence_selection_mutation_0.01/negative/kaks/kaks_sequences.axt.kaks'

# =========================
# -------- FIGURE ---------
# =========================

# =========================
# -------- FIGURE ---------
# =========================

# 1. Use a square figure size
fig = plt.figure(figsize=(12, 12))

# 2. Set height_ratios so all data rows (0,1,2 and 4) are EXACTLY 1.
# The spacer (row 3) is 0.4.
# hspace=0.1 provides a small uniform gap between all rows.
gs = fig.add_gridspec(5, 4, 
                      height_ratios=[1, 1, 1, 0.3, 1], 
                      hspace=0.1, 
                      wspace=0.1)

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

# --- Add the Ghost Plot in Column 3 ---
# This "holds" the 4th column open so plots 0, 1, and 2 don't stretch to fill the row
ghost_ax = fig.add_subplot(gs[4, 3])
ghost_ax.set_visible(False)

# =========================
# -------- PANEL A --------
# =========================

for i in range(3):        # i = row (0, 1, 2)
    for j in range(4):    # j = col (0, 1, 2, 3)
        ax = axes_A[i][j]

        # Load data (using your existing logic)
        fmh_neg = load_fmh(fmh_dnds_negative_files[j][i], 'ref_gene')
        ng_neg = load_ng(ng_dnds_negative_files[i])

        neg, r = combine_and_corr(fmh_neg, ng_neg)

        # Plot the data
        plot_panel(ax, neg, r)
        
        # Formatting Logic:
        # Show ylabel only if leftmost column (j == 0)
        # Show xlabel only if bottom row (i == 2)
        is_left_edge = (j == 0)
        is_bottom_edge = (i == 2)
        
        format_ax(ax, (-0.1, 1.5), (-0.1, 1.5),
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

    # Load data (using your existing logic)
#    fmh_pos = load_fmh(fmh_dnds_positive_files_B[k], 'positive_0.01_ref', upper=8)
    fmh_neg = load_fmh(fmh_dnds_negative_files_B[k], 'negative_0.01_ref', upper=8)

#    ng_pos = load_ng(ng86_dnds_positive_files, upper=8, ref_label='positive_0.01_ref')
    ng_neg = load_ng(ng86_dnds_negative_files, upper=8, ref_label='negative_0.01_ref')

    neg, r = combine_and_corr(fmh_neg, ng_neg)

    # Plot the data
    plot_panel(ax, neg, r)
    
    # Formatting Logic for Row B:
    # ylabel only for the first one (k == 0)
    # xlabel for all of them (since they are the bottom-most row of the figure)
    is_left_edge = (k == 0)
    
    format_ax(ax, (-0.1, 1.5), (-0.1, 1.5),
              xlabel=True, 
              ylabel=is_left_edge)

    # Explicitly hide Y-tick labels for the 2nd and 3rd plots
    if not is_left_edge:
        ax.tick_params(labelleft=False)

# =========================
# -------- FINAL ---------
# =========================

plt.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.08)

fig.text(0.02, 0.95, 'A', fontsize=20, fontweight='bold')
fig.text(0.02, 0.26, 'B', fontsize=20, fontweight='bold')

fig.savefig("/data/jzr5814/repositories/dnds_using_fmh_reproducibles/manuscript_figures/updated_pdf/figure2_sup_negative.pdf",
            bbox_inches='tight')
