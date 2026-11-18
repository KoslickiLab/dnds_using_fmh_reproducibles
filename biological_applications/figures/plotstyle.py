"""Shared look of all figures: fonts, sizes, line widths and colours.

To restyle the figures for a journal, edit the values in this file (figure width, font
family and sizes, colours) and re-run `bash figures/make_all_figures.sh`.
Colours are the colour-blind-safe Okabe-Ito palette."""
import re
import matplotlib

matplotlib.use("Agg")          # draw to files; no display needed
import matplotlib.pyplot as plt

FULL_WIDTH = 7.09              # inches = 180 mm, a full-width journal figure

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],   # the first one installed is used
    "font.size": 7,
    "axes.titlesize": 7.5,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "mathtext.fontset": "custom",       # italic species names use the same sans-serif font
    "mathtext.rm": "sans",
    "mathtext.it": "sans:italic",
    "mathtext.bf": "sans:bold",
    "mathtext.sf": "sans",
    "mathtext.cal": "sans",
    "mathtext.tt": "monospace",
    "mathtext.default": "it",
    "pdf.fonttype": 42,          # keep text as text in the PDF
    "savefig.dpi": 300,
})

# colours with a fixed meaning in every figure
SHARED = "#D55E00"        # near-identical / shared DNA (vermilion)
BACKGROUND = "#8C8C8C"    # everything else (grey)
INK = "#222222"
LIGHT = "#D9D9D9"

ALIGNED_HIGH = SHARED     # alignment bar: at least 99.9% identical
ALIGNED_LOW = "#F0B090"   # alignment bar: 95-99.9% identical

# functional classes of genes (analysis/07_shared_regions.py), in legend order
FUNCTION_COLOURS = {
    "conjugation": "#0072B2",     # blue
    "integration": "#E69F00",     # orange
    "resistance": "#D55E00",      # vermilion
    "defence": "#009E73",         # green
    "transposase": "#CC79A7",     # purple-pink
    "phage": "#56B4E9",           # sky blue
    "replication": "#F0E442",     # yellow
    "other": "#BDBDBD",           # grey
    "hypothetical": "#FFFFFF",    # white
}
FUNCTION_LABELS = {
    "conjugation": "conjugation / mobilization",
    "integration": "integrase / recombinase",
    "resistance": "resistance",
    "defence": "defense systems",
    "transposase": "transposase",
    "phage": "phage",
    "replication": "replication / partition",
    "other": "other function",
    "hypothetical": "hypothetical",
}


def panel_label(ax, letter, dx=-0.02, dy=1.04):
    """Bold panel letter at the top-left corner of an axes."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right")


def italic_species(text):
    """Set species names in italics for matplotlib: 'F. duncaniae A2-165' -> italic 'F. duncaniae' + ' A2-165'.
    Handles abbreviated names anywhere in the text ('other B. intestinihominis strains') and a full name at
    the start ('Roseburia intestinalis ...'; for 'Anaerococcus sp. X' only the genus is italicised)."""
    def italic(words):
        return "$\\it{" + words.replace(" ", "\\ ") + "}$"

    leading = re.match(r"([A-Z][a-z]{3,}) (sp\.|[a-z]{3,})(?= |$)", text)
    if leading:
        genus, species = leading.groups()
        name = genus if species == "sp." else f"{genus} {species}"
        return italic(name) + text[len(name):]
    return re.sub(r"\b([A-Z]\. [a-z]{3,})", lambda match: italic(match.group(1)), text)


def stagger_labels(positions, widths, gap=0.004):
    """Assign each label to a row so that labels in one row do not overlap.

    positions: x of the centre of each label; widths: width of each label (same unit as positions).
    Labels are placed from left to right, each in the lowest row that is free at its position.
    Returns the row number (0 = lowest) of each label, in the order given."""
    rows = [None] * len(positions)
    right_edge_of_row = []
    for index in sorted(range(len(positions)), key=lambda i: positions[i]):
        left, right = positions[index] - widths[index] / 2, positions[index] + widths[index] / 2
        for row, edge in enumerate(right_edge_of_row):
            if left > edge + gap:
                rows[index] = row
                right_edge_of_row[row] = right
                break
        else:
            rows[index] = len(right_edge_of_row)
            right_edge_of_row.append(right)
    return rows


def save(figure, name):
    """Write output/figures/<name>.pdf (for the manuscript) and .png (for a quick look)."""
    import os
    os.makedirs("output/figures", exist_ok=True)
    # no creation date in the files, so that a re-drawn figure is byte-identical to the one in the repository
    figure.savefig(f"output/figures/{name}.pdf", metadata={"CreationDate": None})
    figure.savefig(f"output/figures/{name}.png", dpi=200)
    print(f"wrote output/figures/{name}.pdf and .png")
