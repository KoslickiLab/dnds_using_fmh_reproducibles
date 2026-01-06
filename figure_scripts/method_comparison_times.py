#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --------------------------------------
# Load timing files
# --------------------------------------
path = "/data/jzr5814/repositories/dnds_using_fmh_reproducibles/data/timing_comparison_for_plotting.csv"
df = pd.read_csv(path)

# Explicit sample order
samples = [5, 10, 100, 1000]
x = np.arange(len(samples))

# -----------------------------
# Prepare data for plot 1
# Traditional dN/dS stack plot
# -----------------------------
alignment = df[df["Method"] == "Alignment"].set_index("Sample").loc[samples]["Total (s)"]
ng = df[df["Method"] == "NG Estimate"].set_index("Sample").loc[samples]["Total (s)"]
yn = df[df["Method"] == "YN Estimate"].set_index("Sample").loc[samples]["Total (s)"]

# -----------------------------
# Prepare data for plot 2
# Comparison line plot
# -----------------------------
fmh = df[df["Method"] == "Total FMH"].set_index("Sample").loc[samples]["Total (s)"]
kaks_ng = df[df["Method"] == "Total KaKs NG"].set_index("Sample").loc[samples]["Total (s)"]
kaks_yn = df[df["Method"] == "Total KaKs YN"].set_index("Sample").loc[samples]["Total (s)"]

# -----------------------------
# Plot
# -----------------------------
fig, axes = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(14, 6))
fig.subplots_adjust(wspace=0.05)

# --- Plot 1: Barstack plot ---
axes[0].bar(x, alignment, label="Alignment", color="gray")
axes[0].bar(x, ng, bottom=alignment, label="NG Estimate", color="orange")
axes[0].bar(x, yn, bottom=alignment + ng, label="YN Estimate", color="green")
axes[0].set_title("Traditional dN/dS Walltime")

axes[0].set_xlabel("Sample")
axes[0].set_xticks(x)
axes[0].set_xticklabels(samples)
axes[0].set_ylabel("Time log(seconds)")
axes[0].set_ylim(bottom=1, top=10000)
axes[0].set_yscale("log")
axes[0].legend(loc="upper left")

# --- Plot 2: Line plot ---
offset = 0.05

axes[1].plot(x, fmh, marker="o", label="Total FMH", color="blue")
axes[1].plot(x - offset, kaks_ng, marker="o", label="Total KaKs NG", color="orange")
axes[1].plot(x + offset, kaks_yn, marker="o", label="Total KaKs YN", color="green")

axes[1].set_title("Traditional vs FMH dN/dS Walltimes")
axes[1].set_xlabel("Sample")
axes[1].set_xticks(x)
axes[1].set_xticklabels(samples)
axes[1].set_ylim(bottom=1, top=10000)
axes[1].set_yscale("log")
axes[1].legend()

# --- Save figure ---
output_path = "/data/jzr5814/repositories/dnds_using_fmh_reproducibles/figure_scripts/dnds_runtime_comparison_logscale.pdf"
plt.savefig(output_path, dpi=300, bbox_inches="tight")


#plt.tight_layout()
#plt.show()
