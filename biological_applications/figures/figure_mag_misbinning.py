#!/usr/bin/env python
"""Figure: a contig assigned to the wrong metagenome-assembled genome (MAG).

    usage:  python figures/figure_mag_misbinning.py

Panels
  A  Genes of the Saccharibacteria RGIG2249 scaffold around the region that is shared with
     Methanobrevibacter sp. RGIG2411, and the contig of RGIG2411 that consists of this region.
  B  Read depth along the same stretch in the metagenome from which RGIG2411 was assembled,
     compared with the depth of the archaeal genome and with the depth expected if both
     organisms carried the region (optional steps 16 and 17; results are in the repository).
  C  CheckM contamination estimates of the MAGs and single-cell genomes that carry fully
     shared contigs (pairs of the class "MAG or SAG, foreign contigs").
  D  Number of fully shared contigs in each of these pairs.

Input:  pairs/methanobrevibacter_saccharibacteria/results/ (genes_A.tsv, regions.tsv, region_genes.tsv,
        partner_contigs.tsv, read_depth.tsv, read_depth_windows.tsv),
        all_pairs/results/carrier_genomes.tsv, all_pairs/results/all_pairs_summary.tsv
Output: output/figures/bioapp_mag.pdf and .png
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, Patch
from common import read_settings, pair_accessions
from plotstyle import FULL_WIDTH, SHARED, BACKGROUND, INK, italic_species, save

PAIR = "pairs/methanobrevibacter_saccharibacteria"
FLANK = 23000          # bp shown on either side of the shared region

results = f"{PAIR}/results"
settings = read_settings(PAIR)
name_a, name_b = italic_species(settings["name_a"]), italic_species(settings["name_b"])
gene_labels = dict(item.strip().split("=") for item in settings.get("gene_labels", "").split(";") if "=" in item)
region = pd.read_csv(f"{results}/regions.tsv", sep="\t").iloc[0]
genes = pd.read_csv(f"{results}/genes_A.tsv", sep="\t")
genes = genes[genes.contig == region.contig]
annotated = pd.read_csv(f"{results}/region_genes.tsv", sep="\t")
partner = pd.read_csv(f"{results}/partner_contigs.tsv", sep="\t").iloc[0]
depth = pd.read_csv(f"{results}/read_depth.tsv", sep="\t").set_index("sequence").median_depth
windows = pd.read_csv(f"{results}/read_depth_windows.tsv", sep="\t")
window_start, window_end = max(1, region.start - FLANK), region.end + FLANK

figure = plt.figure(figsize=(FULL_WIDTH, 4.6))
LEFT, WIDTH = 0.085, 0.635            # panels A and B share the x-axis; the space to their right holds the labels of panel B

# ------------------------------------------------------------------ A: genes around the shared region
ax = figure.add_axes([LEFT, 0.775, WIDTH, 0.195])
ax.plot([window_start, window_end], [0, 0], color=INK, linewidth=0.4, zorder=1)
for gene in genes[(genes.end >= window_start) & (genes.start <= window_end)].itertuples():
    start, end = max(gene.start, window_start), min(gene.end, window_end)
    head = min(0.35 * (end - start), 350)
    if gene.strand == "+":
        corners = [(start, -0.2), (end - head, -0.2), (end, 0), (end - head, 0.2), (start, 0.2)]
    else:
        corners = [(end, -0.2), (start + head, -0.2), (start, 0), (start + head, 0.2), (end, 0.2)]
    ax.add_patch(Polygon(corners, facecolor=SHARED if gene.near_identical else "#D0D0D0", edgecolor=INK, linewidth=0.3, zorder=2))
# names of selected genes (settings.ini: gene_labels); names of neighbouring genes are spread apart
labelled = annotated[annotated.locus_tag.isin(gene_labels)].sort_values("start")
gene_x = list((labelled.start + labelled.end) / 2)
SPACING = 0.052 * (window_end - window_start)          # smallest distance between two names
label_x = []
for x in gene_x:
    label_x.append(x if not label_x else max(x, label_x[-1] + SPACING))
for group_start in range(len(label_x)):               # centre each run of shifted names over its genes
    group = [i for i in range(group_start, len(label_x)) if i == group_start or label_x[i] - label_x[i - 1] <= SPACING + 1]
    if group_start == 0 or label_x[group_start] - label_x[group_start - 1] > SPACING + 1:
        run = [group_start]
        while run[-1] + 1 < len(label_x) and label_x[run[-1] + 1] - label_x[run[-1]] <= SPACING + 1:
            run.append(run[-1] + 1)
        shift = (label_x[run[-1]] - gene_x[run[-1]]) / 2
        for i in run:
            label_x[i] -= shift
for gene, x, text_x in zip(labelled.itertuples(), gene_x, label_x):
    ax.plot([x, text_x], [0.2, 0.5], color=INK, linewidth=0.3)
    ax.text(text_x, 0.52, gene_labels[gene.locus_tag], ha="center", va="bottom", fontsize=5.8, style="italic")
# the contig of genome B that consists of the shared region
ax.add_patch(Rectangle((region.start, -0.62), region.end - region.start, 0.14, facecolor=SHARED, linewidth=0, zorder=2))
for x in (region.start, region.end):
    ax.plot([x, x], [-0.72, -0.38], color=INK, linewidth=0.8)
ax.text((region.start + region.end) / 2, -0.78,
        f"contig {partner.contig} of {name_b} ({partner.length_bp / 1000:.1f} kb, {region.identity:.2f}% identical)",
        ha="center", va="top", fontsize=6)
ax.set_xlim(window_start, window_end)
ax.set_ylim(-1.15, 1.1)
ax.axis("off")
ax.legend(handles=[Patch(facecolor=SHARED, edgecolor=INK, linewidth=0.3, label="gene near-identical\nin the archaeal MAG"),
                   Patch(facecolor="#D0D0D0", edgecolor=INK, linewidth=0.3, label="other gene")],
          frameon=False, loc="upper left", bbox_to_anchor=(1.01, 0.75), ncol=1, handlelength=1.1, handletextpad=0.4)
figure.text(0.012, 0.975, "A", fontsize=9, fontweight="bold", va="top")

# ------------------------------------------------------------------ B: read depth along the same stretch
ax = figure.add_axes([LEFT, 0.475, WIDTH, 0.255])
shown = windows[(windows.window_start >= window_start - 1000) & (windows.window_start <= window_end)]
ax.step(shown.window_start, shown.median_depth, where="post", color=INK, linewidth=0.9)
ax.axvspan(region.start, region.end, color=SHARED, alpha=0.12, linewidth=0)
both = depth["expected for the region if both genomes carried it"]
archaeon = depth["genome B (contigs that consist of shared sequence left out)"]
for value, label in ((both, f"expected for the shared region if\nboth organisms carried it ({both:.0f}×)"),
                     (archaeon, f"median depth of the\narchaeal genome ({archaeon:.0f}×)")):
    ax.axhline(value, color=BACKGROUND, linewidth=0.8, linestyle="--")
    ax.text(1.015, value, label, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=6, linespacing=1.15)
ax.set_xlim(window_start, window_end)
ax.set_ylim(0, 1.18 * both)
ax.set_xticks(range(0, int(window_end), 10000), [f"{x // 1000}" for x in range(0, int(window_end), 10000)])
ax.set_xlabel(f"Position on scaffold {region.contig} of {name_a} (kb)")
ax.set_ylabel("Read depth")
figure.text(0.012, 0.745, "B", fontsize=9, fontweight="bold", va="top")

# ------------------------------------------------------------------ C: CheckM contamination of the genomes with wholly shared contigs
ax = figure.add_axes([LEFT, 0.095, 0.36, 0.235])
carriers = pd.read_csv("all_pairs/results/carrier_genomes.tsv", sep="\t")
ax.hist(carriers.checkm_contamination, bins=[0.5 * i for i in range(0, 22)], color=BACKGROUND, linewidth=0)
ax.axvline(5, color=INK, linewidth=0.6, linestyle="--")
ax.text(5.15, ax.get_ylim()[1] * 0.95, "5%", ha="left", va="top", fontsize=6)
_, accession_b = pair_accessions(PAIR)
mine = carriers[carriers.accession == accession_b].checkm_contamination.iloc[0]
ax.annotate("RGIG2411", (mine, 0), xytext=(mine + 0.6, ax.get_ylim()[1] * 0.55), fontsize=6,
            arrowprops=dict(arrowstyle="-", linewidth=0.5, color=INK))
ax.set_xlim(0, 10.5)
ax.set_xlabel("Estimated contamination (CheckM, %)")
ax.set_ylabel("Genomes with\nfully shared contigs")
figure.text(0.012, 0.345, "C", fontsize=9, fontweight="bold", va="top")

# ------------------------------------------------------------------ D: number of fully shared contigs per pair
ax = figure.add_axes([0.60, 0.095, 0.385, 0.235])
pairs = pd.read_csv("all_pairs/results/all_pairs_summary.tsv", sep="\t")
foreign = pairs[pairs["class"] == "MAG or SAG, foreign contigs"]
edges = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]             # bins that double in width
ax.hist(foreign.wholly_shared_contigs, bins=edges, color=BACKGROUND, linewidth=0.5, edgecolor="white")
ax.set_xscale("log")
ax.set_xticks([1, 2, 4, 8, 16, 32, 64, 128, 256], ["1", "2", "4", "8", "16", "32", "64", "128", "256"])
ax.minorticks_off()
ax.set_xlabel("Fully shared contigs per genome pair")
ax.set_ylabel("Genome pairs")
figure.text(0.535, 0.345, "D", fontsize=9, fontweight="bold", va="top")
save(figure, "bioapp_mag")
