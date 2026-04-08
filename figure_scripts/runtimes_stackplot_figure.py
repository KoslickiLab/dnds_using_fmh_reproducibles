import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import matplotlib.patches as mpatches

# Set font style
rcParams['font.family'] = 'serif'

fs=16
rcParams.update({
    'font.size': fs,          # base size
    'axes.labelsize': fs,     # axis labels
    'xtick.labelsize': fs-2,
    'ytick.labelsize': fs-2,
    'legend.fontsize': fs-2,
    'legend.title_fontsize': fs-2
})

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

step_colors = {
    "DNA & Protein Sketches": "#E69F00",  # orange
    "DNA Containment": "#56B4E9",         # sky blue
    "Protein Containment": "#009E73",     # bluish green
    "FracMinHash dN/dS": "#D55E00",       # vermillion
}

# Plot
fig, ax = plt.subplots(figsize=(14, 8))
bar_positions = np.arange(len(pivot_df))

bottom = np.zeros(len(pivot_df))

for step in step_order:
    heights = pivot_df[step].values / 3600  # convert to hours

    ax.bar(
        bar_positions,
        heights,
        bottom=bottom,
        color=step_colors[step],
        edgecolor='black',
        label=step,
        linewidth=0.5
    )

    bottom += heights


# Format x-axis
x_labels = [f'k={k}, t={t}' for k, t in pivot_df.index]
ax.set_xticks(bar_positions)
ax.set_xticklabels(x_labels, rotation=45, ha="right")

# Labels and title
ax.set_ylabel("Runtime (Hours)")


ax.legend(title="Step", loc="upper left", bbox_to_anchor=(1, 1))

plt.tight_layout()

# Save the figure
fig.savefig("/data/jzr5814/repositories/dnds_using_fmh_reproducibles/manuscript_figures/updated_pdf/runtimes_update.png") 
