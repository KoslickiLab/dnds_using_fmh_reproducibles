#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re


# -------------------------------
# 1. Load Ka/Ks Results and Filter
# -------------------------------
kaks_df = pd.read_csv("kaks_results.txt", sep="\t")
# Convert Ka, Ks, and Ka/Ks to numeric (coerce errors to NaN)
kaks_df["Ka"] = pd.to_numeric(kaks_df["Ka"], errors='coerce')
kaks_df["Ks"] = pd.to_numeric(kaks_df["Ks"], errors='coerce')
kaks_df["Ka/Ks"] = pd.to_numeric(kaks_df["Ka/Ks"], errors='coerce')
# Use rows where Ka/Ks > 0.
kaks_df = kaks_df[(kaks_df["Ka"] > 0) & (kaks_df["Ks"] > 0) & (kaks_df["Ka/Ks"] > 0)]

# -------------------------------
# 2. Load Mapping of Query ID to Protein Family
# -------------------------------
mapping_df = pd.read_csv("final_uniprot_mapping.tsv", sep="\t")
protein_family_map = dict(zip(mapping_df["Query ID File 1"].astype(str).str.strip(),
                                mapping_df["Protein Family"].astype(str).str.strip()))

# -------------------------------
# 3. Load the Total Ortholog Counts
# -------------------------------
total_df = pd.read_csv("total_orthologs_between_genomes.tab", sep="\t")
total_map = {}
for idx, row in total_df.iterrows():
    q1 = str(row["Query ID File 1"]).strip()
    q2 = str(row["Query ID File 2"]).strip()
    total_map[(q1, q2)] = row["Total"]

# -------------------------------
# 4. Annotate Ka/Ks Data with Protein Family and Total Count
# -------------------------------
def get_protein_family(seq):
    """Given a sequence field 'ID1_vs_ID2', return the Protein Family for ID1.
       If no mapping exists, return np.nan."""
    parts = seq.split("_vs_")
    id1 = parts[0].strip()
    return protein_family_map.get(id1, np.nan)

def get_total(seq):
    """Given a sequence field 'ID1_vs_ID2', return the 'Total' value from the total_map."""
    parts = seq.split("_vs_")
    if len(parts) < 2:
        return 1
    id1 = parts[0].strip()
    id2 = parts[1].strip()
    return total_map.get((id1, id2), 1)

kaks_df["ProteinFamily"] = kaks_df["Sequence"].apply(get_protein_family)
kaks_df["Total"] = kaks_df["Sequence"].apply(get_total)

# Replace literal "nan" (case-insensitive) with np.nan and remove rows missing a protein family.
kaks_df["ProteinFamily"] = kaks_df["ProteinFamily"].replace(r'(?i)^nan$', np.nan, regex=True)
kaks_df = kaks_df[kaks_df["ProteinFamily"].notna() & (kaks_df["ProteinFamily"].str.strip() != "")]

# -------------------------------
# 5. Prepare Data for Plotting
# -------------------------------
# Calculate the average Ka/Ks for each protein group, sort in descending order
avg_ka_ks = kaks_df.groupby("ProteinFamily")["Ka/Ks"].median().sort_values(ascending=False)
protein_groups = list(avg_ka_ks.index)

def clean_label(label):
    """Remove 'family', 'subfamily', and 'superfamily' (case-insensitive) from the label."""
    return re.sub(r'(?i)(family|subfamily|superfamily)', '', label).strip()

cleaned_groups = [clean_label(x) for x in protein_groups]

method_colors = {
    "NG": "#1f77b4",     # blue
    "LWL": "#ff7f0e",    # orange
    "LPB": "#2ca02c",    # green
    "GY-HKY": "#d62728", # red
    "YN": "#9467bd"      # purple
}

# Define the metrics to plot and their y-axis labels.
metrics = ["Ka", "Ks", "Ka/Ks"]
y_labels = ["dN", "dS", "dN/dS"]

# -------------------------------
# 6. Create Subplots: 3 Rows x 1 Column
# -------------------------------
fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(5, 10))

# Loop over each metric and plot data in its corresponding subplot.
for idx, metric in enumerate(metrics):
    ax = axes[idx]
    for i, pf in enumerate(protein_groups):
        group_data = kaks_df[kaks_df["ProteinFamily"] == pf]
        values = group_data[metric].astype(float)
        
        # Draw a boxplot (without outlier markers).
        bp = ax.boxplot(values, positions=[i], widths=0.6, patch_artist=True,
                        showfliers=False, medianprops=dict(color='black'))
        for box in bp['boxes']:
            box.set_facecolor('white')
            box.set_edgecolor('black')
        
        # Overlay scatter points for each method with a constant marker size.
        for method, color in method_colors.items():
            method_data = group_data[group_data["Method"] == method]
            if method_data.empty:
                continue
            y_vals = method_data[metric].astype(float)
            jitter = np.random.normal(0, 0.08, size=len(method_data))
            x_positions = i + jitter
            ax.scatter(x_positions, y_vals, color=color, alpha=0.9,
                       edgecolors='none', s=20, marker='o', zorder=3)
    
    ax.set_ylabel(y_labels[idx],fontsize=7)
    ax.grid(axis='y', linestyle='--', alpha=0.6)

# Set the x-axis ticks and labels only on the bottom subplot.
axes[-1].set_xticks(range(len(protein_groups)))
axes[-1].set_xticklabels(cleaned_groups, rotation=90, ha='right', fontsize=7)

axes[0].tick_params(axis='y', labelsize=7)
axes[1].tick_params(axis='y', labelsize=7)
axes[2].tick_params(axis='y', labelsize=7)

axes[2].set_xlabel("Orthologs",fontsize=7)

protein_groups = ['rnpA','RNA methyltransferase','RecF','urvA','rpmH','DnaA','LCP','yidC']
ax.set_xticklabels(protein_groups)


plt.tight_layout(rect=[0, 0, 0.85, 1])  # leave room for legend on the right

# -------------------------------
# 7. Add Legend: Method Colors Only
# -------------------------------
from matplotlib.lines import Line2D
method_handles = [Line2D([], [], marker='o', linestyle='None', color=color,
                           label=method, markersize=5)
                  for method, color in method_colors.items()]
# Create custom handles with a specified markersize (e.g., 10)
leg=axes[0].legend(handles=method_handles, title="Method", loc='center left', bbox_to_anchor=(1, 0.75), fontsize=7)
leg.get_title().set_fontsize(7)
plt.rcParams['font.family'] = 'Times New Roman'

plt.savefig("kaks_subplots.pdf")
plt.show()
