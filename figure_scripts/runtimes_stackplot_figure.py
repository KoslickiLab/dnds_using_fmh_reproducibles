import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import matplotlib.patches as mpatches

# Set font style
rcParams['font.family'] = 'serif'

# Load the Excel file
file_path = "../data/log_analysis.xlsx"
xls = pd.ExcelFile(file_path)

# Function to convert hh:mm:ss to seconds
def convert_time(value):
    value = str(value)
    if ":" in value:
        h, m, s = map(int, value.split(":"))
        return h * 3600 + m * 60 + s
    else:
        return int(float(value))

# Read and clean the data
cleaned_data = []
for sheet in xls.sheet_names:
    ksize = int(sheet.split(",")[0].split("=")[1])
    t_threshold = float(sheet.split(",")[1].split("=")[1])
    df = xls.parse(sheet)
    df.columns = ["Step", "Disk Usage (MB)", "Time"]
    df["Time (s)"] = df["Time"].apply(convert_time)
    for _, row in df.iterrows():
        if row["Step"] != "Total":
            cleaned_data.append({
                "Step": row["Step"],
                "Disk Usage (MB)": row["Disk Usage (MB)"],
                "Time (s)": row["Time (s)"],
                "ksize": ksize,
                "t_threshold": t_threshold
            })

df_cleaned = pd.DataFrame(cleaned_data)

# Rename steps and set order
step_mapping = {
    "Sketching": "DNA & Protein Sketches",
    "Pairwise Comparison (DNA)": "DNA Containment",
    "Pairwise Comparison (Protein)": "Protein Containment",
    "dN/dS Estimation": "FracMinHash dN/dS"
}
df_cleaned["Step"] = df_cleaned["Step"].replace(step_mapping)
step_order = ["DNA & Protein Sketches", "DNA Containment", "Protein Containment", "FracMinHash dN/dS"]
df_cleaned["Step"] = pd.Categorical(df_cleaned["Step"], categories=step_order, ordered=True)

# Pivot for stacked bar
pivot_df = df_cleaned.pivot_table(
    index=["ksize", "t_threshold"],
    columns="Step",
    values="Time (s)",
    aggfunc="sum",
    fill_value=0
).sort_index()

# Get unique k-sizes and thresholds
ksizes = sorted(df_cleaned["ksize"].unique())
thresholds = sorted(df_cleaned["t_threshold"].unique())

# Color map and opacity functions
selected_colors = [
    plt.get_cmap("viridis")(0.5),
    plt.get_cmap("viridis")(0.7),
    plt.get_cmap("viridis")(0.9)
]
color_map = lambda i: selected_colors[i % len(selected_colors)]
norm_opacity = plt.Normalize(min(thresholds), max(thresholds))

def adjust_opacity(value, min_opacity=0.3, max_opacity=1.0):
    return max(min_opacity, norm_opacity(value) * (max_opacity - min_opacity) + min_opacity)

# Hatch patterns for each step
hatch_patterns = ["////", "\\\\\\\\", "....", "xxxx"]
step_hatch = {step: hatch_patterns[i % len(hatch_patterns)] for i, step in enumerate(step_order)}

# Plot
fig, ax = plt.subplots(figsize=(14, 8))
bar_positions = np.arange(len(pivot_df))
bottom = np.zeros(len(pivot_df))

for idx, ((ksize, t), row) in enumerate(pivot_df.iterrows()):
    color_idx = ksizes.index(ksize)
    base_color = np.array(color_map(color_idx))
    alpha = adjust_opacity(t)
    rgba_color = (*base_color[:3], alpha)

    for step in step_order:
        height = row[step] / 3600  # Convert to hours
        if height > 0:
            # Draw the bar with the fill color and desired opacity
            ax.bar(
                idx, height, bottom=bottom[idx],
                color=rgba_color, edgecolor='none'
            )
            # Draw the hatch overlay with no fill (ensuring full hatch opacity)
            ax.bar(
                idx, height, bottom=bottom[idx],
                color='none', edgecolor='black', hatch=step_hatch[step]
            )
        bottom[idx] += height


# Format x-axis
x_labels = [f'k={k}, t={t}' for k, t in pivot_df.index]
ax.set_xticks(bar_positions)
ax.set_xticklabels(x_labels, rotation=45, ha="right")

# Labels and title
ax.set_ylabel("Runtime (Hours)")
ax.set_title("FMH dN/dS Runtimes")

# Legend for hatch patterns
legend_handles = [
    mpatches.Patch(facecolor='white', edgecolor='black', hatch=step_hatch[step], label=step)
    for step in step_order
]
ax.legend(handles=legend_handles, title="Step", loc="upper left", bbox_to_anchor=(1, 1))

plt.tight_layout()

# Save the figure
fig.savefig("runtimes.pdf") 
