#!/usr/bin/env python
"""Supplementary figure: how diverged the examined genome pairs are.

    usage:  python figures/figure_supp_pair_divergence.py

One panel per pair of pairs/examined_pairs.tsv: amino-acid identity of the reciprocal-best-hit
orthologs of the two genomes, with the near-identical orthologs in colour.

Input:  <pair folder>/results/orthologs.tsv
Output: output/figures/bioapp_supp_pair_divergence.pdf and .png
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from plotstyle import FULL_WIDTH, SHARED, BACKGROUND, italic_species, save

examined = pd.read_csv("pairs/examined_pairs.tsv", sep="\t")
n_columns = 3
n_rows = -(-len(examined) // n_columns)        # rounded up
height = 1.55 * n_rows + 0.3
figure, axes = plt.subplots(n_rows, n_columns, figsize=(FULL_WIDTH, height), squeeze=False,
                            gridspec_kw=dict(wspace=0.3, hspace=0.95, left=0.07, right=0.985, top=1 - 0.4 / height, bottom=0.42 / height))
for ax in axes.flat[len(examined):]:
    ax.axis("off")

for ax, pair in zip(axes.flat, examined.itertuples()):
    orthologs = pd.read_csv(f"{pair.folder}/results/orthologs.tsv", sep="\t")
    near = orthologs.near_identical
    ax.hist([orthologs[~near].aa_identity, orthologs[near].aa_identity], bins=range(20, 102, 2), stacked=True,
            color=[BACKGROUND, SHARED], linewidth=0)
    ax.set_xlim(20, 100.5)
    first, second = (italic_species(name) for name in pair.name.split(", "))
    ax.set_title(f"{first} and\n{second}", loc="left", fontsize=6.5, linespacing=1.2)
    ax.set_xlabel("Amino-acid identity of orthologs (%)")
for ax in axes[:, 0]:
    ax.set_ylabel("Orthologs")
figure.legend(handles=[Patch(color=SHARED, label="near-identical orthologs"), Patch(color=BACKGROUND, label="other orthologs")],
              frameon=False, loc="lower right", bbox_to_anchor=(0.985, 0.45 / height), ncol=1, handlelength=1.2)
save(figure, "bioapp_supp_pair_divergence")
