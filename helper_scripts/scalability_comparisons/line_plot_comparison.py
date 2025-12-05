#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------
# Load the two timing files
# --------------------------------------
traditional_path = "/data/jzr5814/dnds_scalability_comparison/snakemake_traditional_results/timing_summary_overall.tsv"
modern_path = "/data/jzr5814/dnds_scalability_comparison/snakemake_method_comparison/all_timings_merged.csv"

trad = pd.read_csv(traditional_path, sep="\t")
modern = pd.read_csv(modern_path)

# Ensure pairwise time exists in traditional results
if "time_per_pairwise" not in trad.columns:
    trad["pairwise_count"] = trad["sample"] ** 2
    trad["time_per_pairwise"] = trad["total_seconds"] / trad["pairwise_count"]

# Modern file already contains time_per_pairwise

# Create combined sample labels like "5_rep1"
trad["label"] = trad["sample"].astype(str) + "_rep" + trad["replicate"].astype(str)
modern["label"] = modern["sample"].astype(str) + "_rep" + modern["replicate"].astype(str)

# Sort for consistent ordering
trad = trad.sort_values(["sample", "replicate"])
modern = modern.sort_values(["sample", "replicate"])

# --------------------------------------
# 1. Total runtime comparison (log scale)
# --------------------------------------
plt.figure(figsize=(14, 6))

plt.plot(trad["label"], trad["total_seconds"], marker="o", label="Traditional")
plt.plot(modern["label"], modern["total_time_sec"], marker="o", label="Modern (FMH)")

plt.xticks(rotation=90)
plt.yscale("log")  # LOG TRANSFORM

plt.xlabel("Sample / Replicate")
plt.ylabel("Total Runtime (seconds, log₁₀)")
plt.title("Total Runtime Comparison (Traditional vs FMH)")
plt.legend()
plt.tight_layout()
plt.savefig("runtime_comparison_log.png", dpi=300)
plt.close()

# --------------------------------------
# 2. Time per pairwise comparison (log scale)
# --------------------------------------
plt.figure(figsize=(14, 6))

plt.plot(trad["label"], trad["time_per_pairwise"], marker="o", label="Traditional")
plt.plot(modern["label"], modern["time_per_pairwise"], marker="o", label="Modern (FMH)")

plt.xticks(rotation=90)
plt.yscale("log")  # LOG TRANSFORM

plt.xlabel("Sample / Replicate")
plt.ylabel("Time per Pairwise Comparison (seconds, log₁₀)")
plt.title("Per-Pairwise Runtime Comparison (Traditional vs FMH)")
plt.legend()
plt.tight_layout()
plt.savefig("pairwise_comparison_log.png", dpi=300)
plt.close()

print("Plots saved: runtime_comparison_log.png and pairwise_comparison_log.png")
