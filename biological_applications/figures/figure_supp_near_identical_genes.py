#!/usr/bin/env python
"""Supplementary figure: near-identical genes in the 140 genome pairs of the hierarchical edge
bundling figure.

    usage:  python figures/figure_supp_near_identical_genes.py

Panels
  A  Share of each pair's shared nucleotide 21-mers that lies in near-identical genes, for the
     pairs with a FracMinHash dN/dS of at most 1 and above 1.
  B  FracMinHash dN/dS of the pairs above 1, using all genes and after the near-identical genes
     were removed from both genomes, and of the pairs of related genomes (class "related genomes").
Boxes show the quartiles and the median, whiskers the range of the data within 1.5 times the
interquartile range; every pair is drawn as a point.

Input:  all_pairs/results/all_pairs_summary.tsv
Output: output/figures/bioapp_supp_near_identical.pdf and .png
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotstyle import FULL_WIDTH, SHARED, BACKGROUND, INK, panel_label, save

pairs = pd.read_csv("all_pairs/results/all_pairs_summary.tsv", sep="\t")
above = pairs[pairs.fmh_dnds > 1]
below = pairs[pairs.fmh_dnds <= 1]
related = pairs[pairs["class"] == "related genomes"]
random_numbers = np.random.default_rng(1)          # horizontal jitter of the points


def boxes(ax, groups, colours):
    """One box per group (list of value series), with the individual values as jittered points."""
    values = [group.dropna() for group in groups]
    ax.boxplot(values, positions=range(len(values)), widths=0.5, showfliers=False,
               medianprops=dict(color=INK, linewidth=1.2), boxprops=dict(linewidth=0.7),
               whiskerprops=dict(linewidth=0.7), capprops=dict(linewidth=0.7))
    for position, (group, colour) in enumerate(zip(values, colours)):
        x = position + random_numbers.uniform(-0.17, 0.17, len(group))
        ax.scatter(x, group, s=7, color=colour, alpha=0.7, linewidths=0, zorder=3)


figure, (left, right) = plt.subplots(1, 2, figsize=(FULL_WIDTH, 2.7), gridspec_kw=dict(width_ratios=[2, 3], wspace=0.3,
                                                                                     left=0.08, right=0.985, top=0.92, bottom=0.2))

# A: where the shared k-mers are
boxes(left, [100 * below.near_identical_share, 100 * above.near_identical_share], [BACKGROUND, SHARED])
left.set_xticks([0, 1], [f"at most 1\n(n = {len(below)})", f"above 1\n(n = {len(above)})"])
left.set_xlabel("FracMinHash $d_N/d_S$ of the pair")
left.set_ylabel("Shared 21-mers in\nnear-identical genes (%)")
left.set_ylim(-3, 103)
panel_label(left, "A", dx=-0.2, dy=1.0)

# B: the estimates before and after removal of the near-identical genes
without = above.fmh_dnds_without_near_identical
boxes(right, [above.fmh_dnds, without, related.fmh_dnds], [SHARED, SHARED, BACKGROUND])
right.axhline(1, color=INK, linewidth=0.6, linestyle="--")
right.set_xticks([0, 1, 2], [f"pairs above 1,\nall genes (n = {len(above)})",
                             f"pairs above 1, near-identical\ngenes removed (n = {without.notna().sum()})",
                             f"related genomes,\nall genes (n = {len(related)})"])
right.set_ylabel("FracMinHash $d_N/d_S$")
right.set_ylim(0, 2.1)
panel_label(right, "B", dx=-0.12, dy=1.0)
save(figure, "bioapp_supp_near_identical")
