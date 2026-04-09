#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

# Set font to Times New Roman
rcParams['font.family'] = 'serif'
fs=14
# --------------------------------------
# Load timing files
# --------------------------------------
path = "/data/jzr5814/repositories/dnds_using_fmh_reproducibles/data/timing_comparison_for_plotting.csv"
df = pd.read_csv(path)

# Explicit sample order
samples = [5, 10, 100, 1000]
x = np.arange(len(samples))

fmh = df[df["Method"] == "Total FMH"].set_index("Sample").loc[samples]["Total (s)"]
kaks_ng = df[df["Method"] == "Total KaKs NG"].set_index("Sample").loc[samples]["Total (s)"]
kaks_yn = df[df["Method"] == "Total KaKs YN"].set_index("Sample").loc[samples]["Total (s)"]

fig, ax = plt.subplots(figsize=(7, 6))

offset = 0.05

ax.plot(x, fmh, marker="o", label="FMH dN/dS", color="#448EB9FF")
ax.plot(x - offset, kaks_ng, marker="s", label="NG dN/dS", color="#E69F00")
ax.plot(x + offset, kaks_yn, marker="*", label="YN dN/dS", color="#009E73")

#ax.set_title("Traditional vs FMH dN/dS Walltimes")
ax.set_xlabel("Sample",fontsize=fs)
ax.set_xticks(x)
ax.set_xticklabels(samples,fontsize=fs-2)
ax.set_ylim(bottom=1, top=10000)
ax.set_yscale("log")
ax.tick_params(axis='y', labelsize=fs-2)

ax.legend(loc="upper left")

ax.set_ylabel("Time log(seconds)", fontsize=fs)


# --- Save figure ---
output_path = "/data/jzr5814/repositories/dnds_using_fmh_reproducibles/manuscript_figures/updated_pdf/dnds_runtime_comparison_logscale_line.pdf"
plt.savefig(output_path, dpi=300, bbox_inches="tight")


#plt.tight_layout()
#plt.show()
