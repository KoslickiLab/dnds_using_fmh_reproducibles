"""Shared look of the vignette's figures: the same settings as the main manuscript's figures
(figures/plotstyle.py of the biological applications repository), so the two sets match.
Colours are from the colour-blind-safe Okabe-Ito palette."""
import os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

FULL_WIDTH = 7.09          # inches = 180 mm, a full-width journal figure
HALF_WIDTH = 3.43          # inches = 87 mm

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
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
    "mathtext.fontset": "custom",
    "mathtext.rm": "sans",
    "mathtext.it": "sans:italic",
    "mathtext.bf": "sans:bold",
    "mathtext.sf": "sans",
    "mathtext.default": "it",
    "pdf.fonttype": 42,
    "savefig.dpi": 300,
})

INK = "#222222"
BACKGROUND = "#8C8C8C"
LIGHT = "#D9D9D9"

# one colour per category of isolate pair, in legend order
CATEGORY_COLOURS = {"same sublineage": "#0072B2",                          # blue
                    "same lineage, different sublineage": "#56B4E9",       # sky blue
                    "different lineage": "#D55E00"}                        # vermilion

# one colour per group of genes, in legend order
GROUP_COLOURS = {"essential": "#0072B2",        # blue
                 "all_cds": "#8C8C8C",          # grey: the whole coding genome
                 "nonessential": "#009E73",     # green
                 "pe_ppe": "#E69F00",           # orange
                 "pe_ppe_hq": "#D55E00"}        # vermilion
GROUP_LABELS = {"essential": "essential",
                "all_cds": "whole coding genome",
                "nonessential": "non-essential",
                "pe_ppe": "PE/PPE genes",
                "pe_ppe_hq": "PE/PPE, high-quality\nsegments only"}


def panel_label(ax, letter, dx=-0.02, dy=1.04):
    """Bold panel letter at the top-left corner of an axes."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right")


def box_from_quantiles(ax, position, quantiles, colour, width=0.6, horizontal=False):
    """One box and whisker drawn from (p5, q1, median, q3, p95) instead of from the values."""
    p5, q1, median, q3, p95 = quantiles
    if horizontal:
        ax.add_patch(plt.Rectangle((q1, position - width / 2), q3 - q1, width, facecolor=colour,
                                   edgecolor=INK, linewidth=0.5, zorder=3))
        ax.plot([p5, q1], [position, position], color=INK, linewidth=0.6, zorder=2)
        ax.plot([q3, p95], [position, position], color=INK, linewidth=0.6, zorder=2)
        ax.plot([median, median], [position - width / 2, position + width / 2], color="white",
                linewidth=1.0, zorder=4)
    else:
        ax.add_patch(plt.Rectangle((position - width / 2, q1), width, q3 - q1, facecolor=colour,
                                   edgecolor=INK, linewidth=0.5, zorder=3))
        ax.plot([position, position], [p5, q1], color=INK, linewidth=0.6, zorder=2)
        ax.plot([position, position], [q3, p95], color=INK, linewidth=0.6, zorder=2)
        ax.plot([position - width / 2, position + width / 2], [median, median], color="white",
                linewidth=1.0, zorder=4)


PANEL_SIZE = (3.1, 2.45)   # inches; two panels side by side fit a 6.5-inch text width at their natural size


def _data_points(ax, renderer):
    """Display coordinates of everything plotted in an axes: points along every line (one every two
    pixels), scatter points, the outline and interior of every filled area and box."""
    import numpy as np
    from matplotlib.collections import PathCollection, PolyCollection
    points = []

    def along(vertices):
        for (x0, y0), (x1, y1) in zip(vertices[:-1], vertices[1:]):
            steps = max(2, int(np.hypot(x1 - x0, y1 - y0) / 2))
            points.append(np.column_stack([np.linspace(x0, x1, steps), np.linspace(y0, y1, steps)]))

    for line in ax.get_lines():
        if line.get_visible() and len(line.get_xydata()):
            along(line.get_transform().transform(line.get_xydata()))
    for collection in ax.collections:
        if isinstance(collection, PathCollection):
            points.append(collection.get_offset_transform().transform(collection.get_offsets()))
        elif isinstance(collection, PolyCollection):
            for path in collection.get_paths():
                vertices = collection.get_transform().transform(path.vertices)
                along(vertices)
                x0, y0 = vertices.min(axis=0)
                x1, y1 = vertices.max(axis=0)
                grid = np.mgrid[x0:x1:4, y0:y1:4].reshape(2, -1).T
                if len(grid):
                    points.append(grid[path.transformed(collection.get_transform()).contains_points(grid)])
    for patch in ax.patches:
        box = patch.get_window_extent(renderer)
        grid = np.mgrid[box.x0:box.x1 + 1:2, box.y0:box.y1 + 1:2].reshape(2, -1).T
        points.append(grid)
    points = np.concatenate(points) if points else np.empty((0, 2))
    inside_axes = ax.bbox
    keep = ((points[:, 0] >= inside_axes.x0) & (points[:, 0] <= inside_axes.x1) &
            (points[:, 1] >= inside_axes.y0) & (points[:, 1] <= inside_axes.y1))
    return points[keep]


def overlaps(figure):
    """Every legend or text drawn inside an axes that covers plotted data, as a list of messages."""
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    problems = []
    for ax in figure.axes:
        points = _data_points(ax, renderer)
        boxes = [(f"text '{t.get_text()}'", t.get_window_extent(renderer)) for t in ax.texts]
        if ax.get_legend() is not None:
            boxes.append(("legend", ax.get_legend().get_window_extent(renderer)))
        for name, box in boxes:
            covered = ((points[:, 0] >= box.x0) & (points[:, 0] <= box.x1) &
                       (points[:, 1] >= box.y0) & (points[:, 1] <= box.y1)).sum()
            if covered:
                problems.append(f"{name} covers {covered} plotted points")
    return problems


def save_panel(figure, name):
    """Check that no legend or text covers the data, then write results/figures/panels/<name>.pdf and
    .png (no embedded date, so that reruns match). A panel that fails the check is not written."""
    problems = overlaps(figure)
    os.makedirs("results/figures/panels", exist_ok=True)
    if problems:
        figure.savefig(f"results/figures/panels/{name}_rejected.png", bbox_inches="tight", dpi=150)
        plt.close(figure)
        raise ValueError(f"panel {name}: " + "; ".join(problems) + f" (see results/figures/panels/{name}_rejected.png)")
    if os.path.exists(f"results/figures/panels/{name}_rejected.png"):
        os.remove(f"results/figures/panels/{name}_rejected.png")
    figure.savefig(f"results/figures/panels/{name}.pdf", bbox_inches="tight", metadata={"CreationDate": None})
    figure.savefig(f"results/figures/panels/{name}.png", bbox_inches="tight", dpi=300)
    plt.close(figure)
    print(f"wrote results/figures/panels/{name}.pdf and .png (no overlaps)")


def save(figure, name):
    """Write results/figures/<name>.pdf and .png, without an embedded date so that reruns match."""
    os.makedirs("results/figures", exist_ok=True)
    figure.savefig(f"results/figures/{name}.pdf", bbox_inches="tight",
                   metadata={"CreationDate": None})
    figure.savefig(f"results/figures/{name}.png", bbox_inches="tight", dpi=300)
    plt.close(figure)
    print(f"wrote results/figures/{name}.pdf and .png")
