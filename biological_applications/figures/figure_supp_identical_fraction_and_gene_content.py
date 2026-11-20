#!/usr/bin/env python
"""Supplementary figure: two ways in which a subset of similar sequence raises FracMinHash dN/dS.

    usage:  python figures/figure_supp_identical_fraction_and_gene_content.py

Panels
  A  FracMinHash dN/dS expected when a fraction f of the coding sequence of two genomes is
     identical and the remainder has the containments measured for each examined pair without
     its near-identical genes (lines), with the observed value of each pair (points).
  B  FracMinHash dN/dS of two pairs of complete genomes of one species after a random subset of
     genes was deleted from each genome.

Input:  all_pairs/results/identical_fraction_model.tsv, all_pairs/results/gene_removal.tsv
Output: output/figures/bioapp_supp_identical_fraction.pdf and .png
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from common import fmh_dnds_formula, fmh_dnds_with_identical_fraction
from plotstyle import FULL_WIDTH, INK, panel_label, italic_species, save

PAIR_COLOURS = ["#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#8C8C8C"]

model = pd.read_csv("all_pairs/results/identical_fraction_model.tsv", sep="\t").dropna(subset=["identical_fraction"])
removal = pd.read_csv("all_pairs/results/gene_removal.tsv", sep="\t")

figure, (left, right) = plt.subplots(1, 2, figsize=(FULL_WIDTH, 3.6),
                                     gridspec_kw=dict(wspace=0.3, left=0.075, right=0.985, top=0.94, bottom=0.36))

# A: identical fraction in otherwise diverged genomes
f = np.linspace(0, 0.08, 400)
for pair, colour in zip(model.itertuples(), PAIR_COLOURS):
    left.plot(100 * f, fmh_dnds_with_identical_fraction(f, pair.background_containment_nt, pair.background_containment_aa),
              color=colour, linewidth=1.2, label=" – ".join(italic_species(name) for name in pair.name.split(", ")))
    left.plot(100 * pair.identical_fraction, pair.observed_fmh_dnds, "o", color=colour, markersize=4.5,
              markeredgecolor="white", markeredgewidth=0.5, zorder=3)
left.axhline(1, color=INK, linewidth=0.6, linestyle="--")
left.set_xlim(0, 8)
left.set_ylim(0, 2.6)
left.set_xlabel("Identical fraction of the coding sequence, $f$ (%)")
left.set_ylabel("FracMinHash $d_N/d_S$")
# legend entries show the line and the point of each pair
handles = [Line2D([], [], color=line.get_color(), linewidth=1.2, marker="o", markersize=4, label=line.get_label())
           for line in left.get_lines() if not line.get_label().startswith("_")]
left.legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(-0.14, -0.2), ncol=1, handlelength=1.6, fontsize=6)
panel_label(left, "A", dx=-0.12, dy=1.0)

# B: random gene removal from two complete genomes of one species
for (pair, mine), colour in zip(removal.groupby("pair", sort=False), ["#0072B2", "#E69F00"]):
    means = mine.groupby("kept_percent").fmh_dnds.mean()
    right.scatter(mine.kept_percent, mine.fmh_dnds, s=7, color=colour, alpha=0.55, linewidths=0)
    right.plot(means.index, means.values, color=colour, linewidth=1.2, label=italic_species(pair))
    # expectation: both containments shrink in proportion to the fraction of genes kept
    complete = mine[mine.kept_percent == 100].iloc[0]
    kept = np.linspace(0.5, 1, 100)
    right.plot(100 * kept, fmh_dnds_formula(kept * complete.containment_nt, kept * complete.containment_aa),
               color=colour, linewidth=0.8, linestyle=":")
right.plot([], [], color=INK, linewidth=0.8, linestyle=":", label="expected if both containments shrink\nin proportion to the genes kept")
right.axhline(1, color=INK, linewidth=0.6, linestyle="--")
right.set_yscale("log")
right.set_yticks([0.25, 0.5, 1, 2, 4, 8], ["0.25", "0.5", "1", "2", "4", "8"])
right.minorticks_off()
right.set_ylim(0.2, 10)
right.set_xlim(102, 48)
right.set_xlabel("Genes kept in each genome (%)")
right.set_ylabel("FracMinHash $d_N/d_S$")
right.legend(frameon=False, loc="upper left", bbox_to_anchor=(-0.02, -0.2), ncol=1, handlelength=1.6, fontsize=6)
panel_label(right, "B", dx=-0.12, dy=1.0)
save(figure, "bioapp_supp_identical_fraction")
