import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams


# Set font to Times New Roman
rcParams['font.family'] = 'serif'

# Load the Excel file
file_path = "/data/jzr5814/sourmash_dnds_estimation/thesis_figures/gtdb/log_analysis.xlsx"
xls = pd.ExcelFile(file_path)

# Function to convert time format (hh:mm:ss) to total seconds
def convert_time(value):
    value = str(value)
    if ":" in value:  # Check if it's in hh:mm:ss format
        h, m, s = map(int, value.split(":"))
        return h * 3600 + m * 60 + s
    else:  # Already a numerical value, assume it's in seconds
        return int(float(value))

# Prepare data storage
cleaned_data = []

for sheet in xls.sheet_names:
    ksize = int(sheet.split(",")[0].split("=")[1])  # Extract ksize
    t_threshold = float(sheet.split(",")[1].split("=")[1])  # Extract threshold
    
    # Read and clean sheet data
    df = xls.parse(sheet)
    df.columns = ["Step", "Disk Usage (MB)", "Time"]
    
    # Convert time column properly
    df["Time (s)"] = df["Time"].apply(convert_time)
    
    # Store cleaned data
    for _, row in df.iterrows():
        if row["Step"] not in ["Total"]:  # Exclude total row
            cleaned_data.append({
                "Step": row["Step"],
                "Disk Usage (MB)": row["Disk Usage (MB)"],
                "Time (s)": row["Time (s)"],
                "ksize": ksize,
                "t_threshold": t_threshold
            })

# Convert to DataFrame
df_cleaned = pd.DataFrame(cleaned_data)

# Unique steps and settings
step_mapping = {"Sketching":"DNA & Protein Sketches","dN/dS Estimation": "FracMinHash dN/dS", "Pairwise Comparison (DNA)":"DNA Containment","Pairwise Comparison (Protein)":"Protein Containment"}
df_cleaned["Step"] = df_cleaned["Step"].replace(step_mapping)
df_cleaned["Step"] = pd.Categorical(df_cleaned["Step"], categories=[s for s in df_cleaned["Step"].unique() if s != "FracMinHash dN/dS"] + ["FracMinHash dN/dS"], ordered=True)
df_cleaned = df_cleaned.sort_values("Step")
steps = list(df_cleaned["Step"].unique())
if "FracMinHash dN/dS" in steps:
    steps.append(steps.pop(steps.index("FracMinHash dN/dS")))
ksizes = sorted(df_cleaned["ksize"].unique())
thresholds = sorted(df_cleaned["t_threshold"].unique())

# Define colormap and normalization
cmap = plt.get_cmap("viridis", len(ksizes))
selected_colors = [plt.get_cmap("viridis")(0.0), plt.get_cmap("viridis")(0.4), plt.get_cmap("viridis")(0.7)]
cmap = lambda i: selected_colors[i % len(selected_colors)]

norm_opacity = plt.Normalize(min(thresholds), max(thresholds))

def adjust_opacity(value, min_opacity=0.3, max_opacity=1.0):
    """Adjusts opacity to ensure values are visible."""
    return max(min_opacity, norm_opacity(value) * (max_opacity - min_opacity) + min_opacity)

bar_width = 0.2
spacing = 1.0  # Further increase spacing between categories
indices = np.arange(len(steps)) * (bar_width * len(ksizes) * len(thresholds) + spacing)

# Create first bar plot (Disk Usage)
fig, ax = plt.subplots(figsize=(12, 6))
legend_labels = {}

for i, ksize in enumerate(ksizes):
    for j, t in enumerate(thresholds):
        subset = df_cleaned[(df_cleaned["ksize"] == ksize) & (df_cleaned["t_threshold"] == t)]
        positions = indices + (i * len(thresholds) + j) * bar_width
        label = f"k={ksize}, t={t}"
        opacity = adjust_opacity(t)
        ax.bar(
            positions, np.log10(subset["Disk Usage (MB)"] + 1), width=bar_width,
            color=cmap(i), alpha=opacity, label=label if label not in legend_labels else ""
        )
        legend_labels[label] = True

# Add extended x-ticks to encompass the category
ax.set_xticks(np.concatenate([indices - bar_width, indices + (len(ksizes) * len(thresholds) * bar_width) + bar_width]))
ax.set_xticklabels(["" for _ in range(len(steps) * 2)])
ax.set_xticks(indices + (len(ksizes) * len(thresholds) * bar_width) / 2, minor=True)
ax.set_xticklabels(steps, minor=True, rotation=360, ha="center" )
ax.tick_params(axis='x', which='major', length=20)
ax.set_ylabel("Disk Usage Log10(MB)")
ax.set_title("FMH dN/dS Disk Usage")
ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()


# Save the figure
fig.savefig("disk_usage.pdf") 