"""The vignette's figures, one file per panel, from the FracMinHash dN/dS tool's output only.

    usage:  python figures/42_figures.py

Each panel is drawn on its own and checked: a panel in which a legend or a text covers any plotted
point, line or box is not written (plotstyle.save_panel). The panels are put together, with their
letters, in tex/mtb_figures.tex. Titles and sample sizes go in the captions, not in the panels.

Reads results/tool/ (step 40), results/trends_tool/ (steps 41, 44), results/conventional_tool/
(step 33) and results/smoke/ (step 12). Writes results/figures/panels/:
    mtb_main_a   whole-coding-genome values by how closely the two isolates are related
    mtb_main_b   the groups of genes, with the per-gene values of Gan et al.
    mtb_main_c   FracMinHash dN/dS against seqinr on the same sequences
    mtb_main_d   whole-coding-genome values against the divergence of the two isolates
    mtb_scale_a  nucleotide against protein containment: draft assemblies and anchored sequences
    mtb_scale_b  agreement of sketched values with those from every k-mer, by scale factor
    mtb_scale_c  median value of each group, by scale factor
    mtb_scale_d  how often the orderings of Gan et al. hold, by scale factor
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plotstyle import (CATEGORY_COLOURS, GROUP_COLOURS, GROUP_LABELS, INK, LIGHT, PANEL_SIZE,
                       box_from_quantiles, save_panel)
import matplotlib.pyplot as plt

MTB = os.getcwd()
GROUP_ORDER = ["essential", "all_cds", "nonessential", "pe_ppe", "pe_ppe_hq"]
# per-gene pN/pS reported by Gan et al. 2025; their PE/PPE value is over the high-quality segments
PUBLISHED = {"essential": 0.69, "nonessential": 0.84, "pe_ppe_hq": 0.90}
TRENDS = f"{MTB}/results/trends_tool"
SHORT_LABELS = {"essential": "essential", "all_cds": "whole coding genome", "nonessential": "non-essential",
                "pe_ppe": "PE/PPE genes", "pe_ppe_hq": "PE/PPE segments"}
COLUMNS = ["same_sublineage", "same_lineage", "different_lineage"]
SCALE_TICKS = [1, 2, 3, 5, 10, 100, 1000]
SEQINR_LIMIT = 1.2               # upper limit of both axes of main panel c
DIVERGENCE_MIN_PAIRS = 100_000   # divergence bins with fewer pairs are not drawn in main panel d
# what the captions say about the panels, written to results/figures/parameters.tsv at the end
PARAMETERS = []


def read(path):
    return pd.read_csv(path, sep="\t")


def histogram(group, scale=1):
    return read(f"{MTB}/results/tool/{group}/scaled{scale}/histogram.tsv")


def scale_axis(ax):
    ax.set_xscale("log")
    ax.set_xticks(SCALE_TICKS)
    ax.set_xticklabels([str(s) for s in SCALE_TICKS])
    ax.minorticks_off()
    ax.set_xlim(0.8, 1300)
    ax.set_xlabel("scale factor")


def panel():
    return plt.subplots(figsize=PANEL_SIZE)


# ------------------------------------------------------------------ Figure 1
# (a) the whole coding genome, by how closely the two isolates are related
figure, ax = panel()
table = histogram("all_cds").iloc[1:-1]
centres = 0.5 * (table.bin_low + table.bin_high)
PARAMETERS += [("histogram_bin_width", round(float((table.bin_high - table.bin_low).iloc[0]), 6),
                "width of the bins of main panel a")]
labels = {"same_sublineage": ("same sublineage", "same sublineage"),
          "same_lineage": ("same lineage, different sublineage", "same lineage, different sublineages"),
          "different_lineage": ("different lineage", "different lineages")}
highest = 0
for column, (category, label) in labels.items():
    share = table[column] / table[column].sum()
    highest = max(highest, share.max())
    ax.plot(centres, share, color=CATEGORY_COLOURS[category], linewidth=1.0, label=label)
ax.axvline(1.0, color=INK, linewidth=0.6, linestyle=(0, (3, 2)))
ax.set_xlim(0.2, 1.05)
ax.set_ylim(0, highest * 1.55)
ax.set_xlabel("FracMinHash dN/dS, whole coding genome")
ax.set_ylabel("share of pairs of isolates")
ax.legend(frameon=False, handlelength=1.4, loc="upper left")
save_panel(figure, "mtb_main_a")

# (b) the groups of genes (narrower plot, to leave room for the names of the groups)
figure, ax = plt.subplots(figsize=(PANEL_SIZE[0] - 0.6, PANEL_SIZE[1]))
rows = [g for g in GROUP_ORDER]
groups = read(f"{TRENDS}/groups.tsv").set_index("group")      # exact quantiles over every pair (step 40b)
for index, group in enumerate(rows):
    box = groups.loc[group, ["p5", "q1", "median", "q3", "p95"]].to_numpy(dtype=float)
    box_from_quantiles(ax, -index, box, GROUP_COLOURS[group], width=0.55, horizontal=True)
    if group in PUBLISHED:
        ax.plot([PUBLISHED[group]], [-index], marker="D", markersize=3.4, markerfacecolor="white",
                markeredgecolor=INK, markeredgewidth=0.7, linestyle="none", zorder=5,
                label="Gan et al. 2025, per-gene pN/pS" if group == "essential" else None)
ax.vlines(1.0, -len(rows) + 0.5, 0.5, color=INK, linewidth=0.6, linestyle=(0, (3, 2)))   # across the rows only
ax.set_yticks([-i for i in range(len(rows))])
ax.set_yticklabels([GROUP_LABELS[g] for g in rows])
ax.set_ylim(-len(rows) + 0.4, 1.3)
ax.set_xlabel("FracMinHash dN/dS")
ax.tick_params(axis="y", length=0)
ax.legend(frameon=False, handlelength=1.0, loc="upper left")
save_panel(figure, "mtb_main_b")

# (c) against seqinr, same sequences. Grey points: the whole coding genome, one per pair. Crosses:
# for each group of genes, the median of both estimates with their quartiles (colours as in b).
# A few pairs with almost no synonymous differences reach values far above the rest; the axes stop at
# SEQINR_LIMIT and the number of whole-genome pairs beyond it is written out for the caption.
figure, ax = panel()
whole = read(f"{MTB}/results/conventional_tool/all_cds.tsv").dropna(subset=["seqinr_dnds", "fmh_dnds"])
beyond = (whole.seqinr_dnds > SEQINR_LIMIT) | (whole.fmh_dnds > SEQINR_LIMIT) | (whole.fmh_dnds < 0)
ax.scatter(whole.seqinr_dnds[~beyond], whole.fmh_dnds[~beyond], s=1.0, linewidths=0, alpha=0.3,
           color=GROUP_COLOURS["all_cds"], rasterized=True)
for group in GROUP_ORDER:
    pairs = read(f"{MTB}/results/conventional_tool/{group}.tsv").dropna(subset=["seqinr_dnds", "fmh_dnds"])
    x1, x, x3 = pairs.seqinr_dnds.quantile([0.25, 0.5, 0.75])
    y1, y, y3 = pairs.fmh_dnds.quantile([0.25, 0.5, 0.75])
    ax.errorbar(x, y, xerr=[[x - x1], [x3 - x]], yerr=[[y - y1], [y3 - y]], fmt="o", markersize=3.6,
                color=GROUP_COLOURS[group], markeredgecolor=INK, markeredgewidth=0.5, elinewidth=1.0,
                capsize=0, zorder=5)
ax.plot([0, SEQINR_LIMIT], [0, SEQINR_LIMIT], color=INK, linewidth=0.6, linestyle=(0, (3, 2)))
ax.set_xlim(0, SEQINR_LIMIT)
ax.set_ylim(0, SEQINR_LIMIT)
ax.set_aspect("equal")
ax.set_xlabel("seqinr dN/dS")
ax.set_ylabel("FracMinHash dN/dS")
save_panel(figure, "mtb_main_c")
PARAMETERS += [("seqinr_limit", SEQINR_LIMIT, "upper limit of both axes of main panel c"),
               ("seqinr_pairs_outside", int(beyond.sum()), "whole-genome pairs beyond those axes, not drawn")]

# (d) against the divergence of the two isolates of a pair
figure, ax = panel()
divergence = read(f"{TRENDS}/divergence.tsv")
divergence = divergence[divergence.pairs >= DIVERGENCE_MIN_PAIRS]
PARAMETERS += [("divergence_min_pairs", DIVERGENCE_MIN_PAIRS, "fewest pairs in a divergence bin drawn in main panel d")]
middle = 0.5 * (divergence.divergence_low + divergence.divergence_high)
ax.fill_between(middle * 100, divergence.q1, divergence.q3, color=LIGHT, linewidth=0)
ax.plot(middle * 100, divergence["median"], color="#0072B2", linewidth=1.2)
ax.set_xlabel("divergence between the isolates of a pair, $1 - C_{nt}^{1/21}$ (%)")
ax.set_ylabel("FracMinHash dN/dS, whole coding genome")
save_panel(figure, "mtb_main_d")

# ------------------------------------------------------------------ Figure 2
scales = read(f"{TRENDS}/scales.tsv")
ordering = read(f"{TRENDS}/ordering.tsv")

# (a) nucleotide against protein containment, two ways
figure, ax = panel()
smoke = read(f"{MTB}/results/smoke/pairs.tsv")
ax.scatter(1 - smoke.C_nt_1, 1 - smoke.C_aa_1, s=2.0, linewidths=0, alpha=0.5, color="#D55E00",
           label="draft assemblies,\nown gene calls", rasterized=True)
anchored = read(f"{TRENDS}/containment_sample.tsv")         # pairs drawn at random in step 41
PARAMETERS += [("anchored_points", len(anchored), "anchored pairs drawn in scale panel a")]
ax.scatter(1 - anchored.C_nt, 1 - anchored.C_aa, s=2.0, linewidths=0, alpha=0.5, color="#0072B2",
           label="anchored to one reference", rasterized=True)
top = 1.05 * max((1 - smoke.C_nt_1).max(), (1 - smoke.C_aa_1).max())
ax.plot([0, top], [0, top], color=INK, linewidth=0.6, linestyle=(0, (3, 2)))
ax.set_xlim(0, top)
ax.set_ylim(0, top)
ax.set_xlabel("$1 - C_{nt}$")
ax.set_ylabel("$1 - C_{aa}$")
ax.legend(frameon=False, handlelength=1.0, markerscale=3, loc="lower right")
save_panel(figure, "mtb_scale_a")

# (b) agreement with every k-mer, (c) median, by scale factor
for name, column, ylabel in (("mtb_scale_b", "spearman_with_scale_1", "rank correlation with every $k$-mer"),
                             ("mtb_scale_c", "median", "median FracMinHash dN/dS")):
    figure, ax = panel()
    for group in GROUP_ORDER:
        rows = scales[scales.group == group].sort_values("scale")
        ax.plot(rows.scale, rows[column], marker="o", markersize=2.6, linewidth=1.0, color=GROUP_COLOURS[group],
                label=SHORT_LABELS[group])
    scale_axis(ax)
    ax.set_ylabel(ylabel)
    if name == "mtb_scale_b":
        ax.set_ylim(-0.4, 1.05)
        ax.legend(frameon=False, handlelength=1.4, loc="lower left", fontsize=6)
    else:
        ax.axhline(1.0, color=INK, linewidth=0.6, linestyle=(0, (3, 2)))
        ax.set_ylim(-0.05, 1.1)
    save_panel(figure, name)

# (d) the orderings of Gan et al., pair by pair
figure, ax = panel()
styles = {("essential", "nonessential"): ("essential < non-essential", "#009E73"),
          ("essential", "pe_ppe_hq"): ("essential < PE/PPE segments", "#0072B2"),
          ("nonessential", "pe_ppe_hq"): ("non-essential < PE/PPE segments", "#D55E00")}
for (lower, higher), (label, colour) in styles.items():
    rows = ordering[(ordering.lower == lower) & (ordering.higher == higher)].sort_values("scale")
    ax.plot(rows.scale, rows.share_holding, marker="o", markersize=2.6, linewidth=1.0, color=colour, label=label)
ax.axhline(0.5, color=INK, linewidth=0.6, linestyle=(0, (3, 2)))
scale_axis(ax)
ax.set_ylim(0.0, 1.32)                   # room above the data for the legend
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_ylabel("share of pairs in which it holds")
ax.legend(frameon=False, handlelength=1.4, loc="upper right", fontsize=6)
save_panel(figure, "mtb_scale_d")

os.makedirs(f"{MTB}/results/figures", exist_ok=True)
pd.DataFrame(PARAMETERS, columns=["name", "value", "meaning"]).to_csv(
    f"{MTB}/results/figures/parameters.tsv", sep="\t", index=False)
print("wrote results/figures/parameters.tsv")
