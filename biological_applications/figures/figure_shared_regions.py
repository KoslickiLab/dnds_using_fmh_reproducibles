#!/usr/bin/env python
"""Figure of one genome pair: where the shared DNA lies, what it encodes, and which other
genomes carry it.

    usage:  python figures/figure_shared_regions.py pairs/<name>

Panels
  A  For every gene of genome A, the fraction of its 21-mers that also occur in genome B, along
     genome A. The shared regions (E1, E2, ...) are named above the plot.
  B  Gene maps of selected regions (settings.ini: regions_to_draw), coloured by function, with a
     bar underneath showing what aligns to genome B and at what identity.
  C  Presence of every shared region in genome B and in the other genomes of genomes.tsv.

Input:  the result tables in pairs/<name>/results/ (no genome data needed)
Output: output/figures/<figure_name of settings.ini>.pdf and .png

The look of the figure (fonts, colours, width) is set in figures/plotstyle.py; the sizes of the
panels are set by the constants below.
"""
import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, Patch
from matplotlib.lines import Line2D
from common import read_genomes, read_settings
from plotstyle import (FULL_WIDTH, SHARED, BACKGROUND, INK, LIGHT, ALIGNED_HIGH, ALIGNED_LOW, FUNCTION_COLOURS,
                       FUNCTION_LABELS, panel_label, italic_species, stagger_labels, save)

# ---- sizes, in inches
TRACK_HEIGHT = 1.05           # panel A, without the region names above it
LABEL_ROW_HEIGHT = 0.13       # one row of region names above panel A
GENE_MAP_ROW_HEIGHT = 0.43    # one region in panel B
LEGEND_HEIGHT = 0.50          # legends under panel B
HEATMAP_ROW_HEIGHT = 0.15     # one genome in panel C
HEATMAP_CELL_WIDTH = 0.30     # one region in panel C (made narrower when there are many regions)
LEFT_MARGIN = 0.60            # room for the y-axis label of panel A
GROUP_LABEL_WIDTH = 1.50      # room for the group names on the right of panel C

pair_dir = sys.argv[1]
results = f"{pair_dir}/results"
settings = read_settings(pair_dir)
genomes = read_genomes(pair_dir)
name_a, name_b = italic_species(settings["name_a"]), italic_species(settings["name_b"])

genes_a = pd.read_csv(f"{results}/genes_A.tsv", sep="\t")
contigs = pd.read_csv(f"{results}/contigs_A.tsv", sep="\t")
regions = pd.read_csv(f"{results}/regions.tsv", sep="\t")
region_genes = pd.read_csv(f"{results}/region_genes.tsv", sep="\t")
blocks = pd.read_csv(f"{results}/blocks.tsv", sep="\t")
presence = pd.read_csv(f"{results}/presence.tsv", sep="\t")

regions_to_draw = [r.strip() for r in settings["regions_to_draw"].split(",")]
# optional names above single genes in panel B, written in settings.ini as  locus_tag=label; locus_tag=label
gene_labels = dict(item.strip().split("=") for item in settings.get("gene_labels", "").split(";") if "=" in item)
others = genomes[genomes.role != "A"].reset_index(drop=True)      # rows of panel C

# ---- positions along genome A: contigs are laid end to end, longest first
contigs["offset"] = contigs.length_bp.cumsum() - contigs.length_bp
offset = dict(zip(contigs.contig, contigs.offset))
genome_length = contigs.length_bp.sum()
plot_width = FULL_WIDTH - LEFT_MARGIN - 0.1                         # width of panels A and B, inches

# ---- names of the regions above panel A, in as many rows as needed to keep them apart
region_x = [(offset[r.contig] + (r.start + r.end) / 2) / genome_length for r in regions.itertuples()]      # 0..1
label_width = [(0.055 + 0.045 * len(name)) / plot_width for name in regions.region]                           # 0..1
label_row = stagger_labels([min(max(x, w / 2), 1 - w / 2) for x, w in zip(region_x, label_width)], label_width)
n_label_rows = max(label_row) + 1

# ---- layout: one column of panels; every height in inches
heights = dict(track=TRACK_HEIGHT + LABEL_ROW_HEIGHT * n_label_rows + 0.1,
               track_axis=0.62,                                     # x-axis label and legend of panel A
               maps=GENE_MAP_ROW_HEIGHT * len(regions_to_draw) + 0.22,
               legends=LEGEND_HEIGHT,
               heatmap_header=0.42,                                 # region names and sizes above panel C
               heatmap=HEATMAP_ROW_HEIGHT * len(others),
               colourbar=0.45)
figure_height = sum(heights.values()) + 0.1
figure = plt.figure(figsize=(FULL_WIDTH, figure_height))


def axes_at(top, height, left=LEFT_MARGIN, width=plot_width):
    """Axes placed by inches from the top-left corner of the figure."""
    return figure.add_axes([left / FULL_WIDTH, 1 - (top + height) / figure_height, width / FULL_WIDTH, height / figure_height])


# ------------------------------------------------------------------ A: shared k-mers along genome A
top = 0.1
ax = axes_at(top, heights["track"])
position = (genes_a.contig.map(offset) + (genes_a.start + genes_a.end) / 2) / 1e6
near = genes_a.near_identical
ax.scatter(position[~near], genes_a.fraction_shared[~near], s=2.5, color=BACKGROUND, linewidths=0, rasterized=True)
ax.scatter(position[near], genes_a.fraction_shared[near], s=4, color=SHARED, linewidths=0, rasterized=True)
if len(contigs) > 1:                                  # mark contig boundaries of a draft genome
    for contig in contigs.itertuples():
        if contig.offset > 0 and contig.length_bp > 0.01 * genome_length:
            ax.plot([contig.offset / 1e6] * 2, [-0.03, 1.03], color=LIGHT, linewidth=0.4, zorder=0)
# region names, each joined to its region by a thin line
units_per_inch = 1.06 / TRACK_HEIGHT                  # y-axis units per inch
for region, x, width, row in zip(regions.itertuples(), region_x, label_width, label_row):
    y = 1.06 + (0.05 + LABEL_ROW_HEIGHT * row) * units_per_inch
    label_x = min(max(x, width / 2), 1 - width / 2)          # keep names at the two ends of the genome inside the panel
    ax.plot([x * genome_length / 1e6, label_x * genome_length / 1e6], [1.04, y], color=BACKGROUND, linewidth=0.4)
    ax.text(label_x * genome_length / 1e6, y, region.region, ha="center", va="bottom", fontsize=6,
            fontweight="bold" if region.region in regions_to_draw else "normal")
ax.set_xlim(0, genome_length / 1e6)
ax.set_ylim(-0.03, 1.06 + (0.05 + LABEL_ROW_HEIGHT * n_label_rows + 0.05) * units_per_inch)
ax.set_yticks([0, 0.5, 1])
ax.spines["left"].set_bounds(0, 1)
ax.set_xlabel(f"Position in the {name_a} genome (Mb)" + ("" if len(contigs) == 1 else "; contigs joined end to end, longest first"),
              loc="left")
ax.set_ylabel("Fraction of the gene's\n21-mers shared", y=0.5 / ax.get_ylim()[1], ha="center")
ax.legend(handles=[Line2D([], [], marker="o", linestyle="", color=SHARED, markersize=3.5,
                          label=f"near-identical gene (n = {int(near.sum())})"),
                   Line2D([], [], marker="o", linestyle="", color=BACKGROUND, markersize=3,
                          label=f"other gene (n = {int((~near).sum()):,})")],
          frameon=False, loc="upper right", bbox_to_anchor=(1.0, -0.36 / heights["track"]), ncol=2,
          handletextpad=0.2, columnspacing=1.2, borderaxespad=0)
figure.text(0.012, 1 - (top + 0.02) / figure_height, "A", fontsize=9, fontweight="bold", va="top")

# ------------------------------------------------------------------ B: gene maps of selected regions
top += heights["track"] + heights["track_axis"]
ax = axes_at(top, heights["maps"])
drawn = regions.set_index("region").loc[regions_to_draw]
longest = drawn.span_bp.max()
for row, (name, region) in enumerate(drawn.iterrows()):
    y = -row
    ax.plot([0, region.span_bp], [y, y], color=INK, linewidth=0.4, zorder=1)
    last_label_x = None          # position of the previous gene label in this row, to stagger close neighbours
    label_level = 0
    for gene in region_genes[region_genes.region == name].sort_values("start").itertuples():
        start = max(gene.start, region.start) - region.start
        end = min(gene.end, region.end) - region.start
        head = min(0.35 * (end - start), 0.006 * longest)       # length of the arrow head
        if gene.strand == "+":
            corners = [(start, y - 0.2), (end - head, y - 0.2), (end, y), (end - head, y + 0.2), (start, y + 0.2)]
        else:
            corners = [(end, y - 0.2), (start + head, y - 0.2), (start, y), (start + head, y + 0.2), (end, y + 0.2)]
        ax.add_patch(Polygon(corners, facecolor=FUNCTION_COLOURS[gene.function], edgecolor=INK, linewidth=0.3, zorder=2))
        if gene.locus_tag in gene_labels:
            x = (start + end) / 2
            close_to_previous = last_label_x is not None and x - last_label_x < 0.12 * longest
            label_level = 1 - label_level if close_to_previous else 0
            ax.plot([x, x], [y + 0.2, y + 0.26 + 0.2 * label_level], color=INK, linewidth=0.3)
            ax.text(x, y + 0.27 + 0.2 * label_level, gene_labels[gene.locus_tag], ha="center", va="bottom", fontsize=5.5)
            last_label_x = x
    # bar under the genes: what aligns to genome B (dark = at least 99.9% identical, light = 95-99.9%)
    for block in blocks[blocks.region == name].itertuples():
        ax.add_patch(Rectangle((block.start_ref - region.start, y - 0.37), block.end_ref - block.start_ref, 0.09,
                               facecolor=ALIGNED_HIGH if block.identity >= 99.9 else ALIGNED_LOW, linewidth=0, zorder=2))
    ax.text(-0.012 * longest, y, name, ha="right", va="center", fontweight="bold")
    ax.text(region.span_bp + 0.008 * longest, y,
            f"{region.span_bp / 1000:.0f} kb, {min(region.identity, 99.99):.2f}%", ha="left", va="center", fontsize=6)
# scale bar
scale_y = -len(drawn) + 0.42
ax.plot([0, 10000], [scale_y, scale_y], color=INK, linewidth=1)
ax.text(10000 + 0.006 * longest, scale_y, "10 kb", ha="left", va="center", fontsize=6)
ax.set_xlim(-0.03 * longest, 1.13 * longest)
ax.set_ylim(-len(drawn) + 0.3, 0.6)
ax.axis("off")
figure.text(0.012, 1 - (top - 0.02) / figure_height, "B", fontsize=9, fontweight="bold", va="top")

# legends of panel B, in their own strip
top += heights["maps"]
ax = axes_at(top, heights["legends"])
ax.axis("off")
present_functions = [f for f in FUNCTION_COLOURS if (region_genes[region_genes.region.isin(regions_to_draw)].function == f).any()]
handles = [Patch(facecolor=FUNCTION_COLOURS[f], edgecolor=INK, linewidth=0.3, label=FUNCTION_LABELS[f]) for f in present_functions]
gene_legend = ax.legend(handles=handles, frameon=False, ncol=3, loc="upper left", bbox_to_anchor=(0.0, 1.0),
                        handlelength=1.1, columnspacing=1.2, handletextpad=0.4, borderaxespad=0)
ax.add_artist(gene_legend)
ax.legend(handles=[Patch(facecolor=ALIGNED_HIGH, label="≥99.9% identical"), Patch(facecolor=ALIGNED_LOW, label="95–99.9% identical")],
          frameon=False, loc="upper right", bbox_to_anchor=(1.0, 1.0), handlelength=1.1, handletextpad=0.4, borderaxespad=0,
          title=f"Aligned to {name_b}", alignment="left")

# ------------------------------------------------------------------ C: presence in other genomes
top += heights["legends"] + heights["heatmap_header"]
region_names = list(regions.region)
heatmap_left = 0.55 + 0.056 * others.label.str.len().max()        # room for the genome names (about 0.056 inch per letter)
cell_width = min(HEATMAP_CELL_WIDTH, (FULL_WIDTH - heatmap_left - GROUP_LABEL_WIDTH) / len(region_names))
ax = axes_at(top, heights["heatmap"], left=heatmap_left, width=cell_width * len(region_names))
matrix = presence.pivot(index="accession", columns="region", values="covered").loc[others.accession, region_names]
image = ax.imshow(matrix.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(len(region_names)), [f"{r}\n{size:.0f}" for r, size in zip(region_names, regions.span_bp / 1000)],
              fontsize=5.5 if cell_width < 0.25 else 6.5)
ax.text(-0.01, 1.0, "region\nsize (kb)", transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5, linespacing=1.25)
ax.set_yticks(range(len(others)), [italic_species(label) for label in others.label])
ax.tick_params(length=0)
ax.xaxis.tick_top()
for spine in ax.spines.values():
    spine.set_visible(False)
ax.set_xticks([x - 0.5 for x in range(1, len(region_names))], minor=True)       # thin white gaps between cells
ax.set_yticks([y - 0.5 for y in range(1, len(others))], minor=True)
ax.grid(which="minor", color="white", linewidth=1.0)
ax.tick_params(which="minor", length=0)
# a line and a label on the right for each group of genomes
for group, rows in others.groupby("group", sort=False):
    if group == "pair":
        group = "other genome of the pair"
    first, last = rows.index.min(), rows.index.max()
    x = len(region_names) - 0.5 + 0.06 / cell_width
    ax.plot([x, x], [first - 0.4, last + 0.4], color=INK, linewidth=0.7, clip_on=False)
    ax.text(x + 0.05 / cell_width, (first + last) / 2, italic_species(group), va="center", ha="left", fontsize=6, clip_on=False)
colourbar_axes = axes_at(top + heights["heatmap"] + 0.12, 0.07, left=heatmap_left, width=0.9)
colourbar = figure.colorbar(image, cax=colourbar_axes, orientation="horizontal", ticks=[0, 0.5, 1])
colourbar.outline.set_linewidth(0.4)
colourbar.ax.tick_params(labelsize=6, length=2)
colourbar_axes.text(1.08, 0.5, "fraction of the region present (≥90% identity)", transform=colourbar_axes.transAxes,
                    va="center", ha="left", fontsize=6)
figure.text(0.012, 1 - (top - heights["heatmap_header"] + 0.02) / figure_height, "C", fontsize=9, fontweight="bold", va="top")

save(figure, settings["figure_name"])
