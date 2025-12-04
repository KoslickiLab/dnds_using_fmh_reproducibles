#!/usr/bin/env python3
import os
import argparse
import pandas as pd

def normalize_sample(df, sample_col):
    # Rename original column to preserve full identifier
    df = df.rename(columns={sample_col: f"{sample_col}_replicate"})

    # Extract the sample number (e.g., sample5_rep1 → 5)
    df["sample"] = (
        df[f"{sample_col}_replicate"]
        .str.extract(r"sample(\d+)", expand=False)
        .astype(int)
    )

    # Extract the replicate number (e.g., sample5_rep1 → 1)
    df["replicate"] = (
        df[f"{sample_col}_replicate"]
        .str.extract(r"rep(\d+)", expand=False)
        .astype(int)
    )

    return df




def main():
    parser = argparse.ArgumentParser(description="Merge all timing logs into one dataset.")
    parser.add_argument("--base", required=True, help="Base directory containing timing logs.")
    args = parser.parse_args()

    base = args.base

    # --- File paths ---
    sig_path = os.path.join(base, "signature_generation_times.csv")
    cont_path = os.path.join(base, "containment_times.csv")
    fmh_path = os.path.join(base, "fmh_omega_times.csv")
    out_path = os.path.join(base, "all_timings_merged.csv")

    # --- Load ---
    sig_df = pd.read_csv(sig_path)
    cont_df = pd.read_csv(cont_path)
    fmh_df = pd.read_csv(fmh_path)

    #debugging normalize sample
#    print("DEBUG sample column dtype:", sig_df["sample"].dtype)
#    print("DEBUG sample column head:")
#    print(sig_df["sample"].head())

    # --- Normalize sample identifiers ---
    sig_df = normalize_sample(sig_df, "sample")
    # For sketches: rename sample5_rep1 column to something else
    if "sample" in sig_df.columns and sig_df["sample"].dtype == object:
        sig_df.rename(columns={"sample": "sample_raw"}, inplace=True)
        sig_df = normalize_sample(sig_df, "sample_raw")

    # Containment DF already has sample + replicate
    cont_df["sample"] = cont_df["sample"].astype(int)
    cont_df["replicate"] = cont_df["replicate"].astype(int)

    # FMH DF already has sample + replicate
    fmh_df["sample"] = fmh_df["sample"].astype(int)
    fmh_df["replicate"] = fmh_df["replicate"].astype(int)

    # --- Pivot containment so DNA & protein times become columns ---
    cont_pivot = cont_df.pivot_table(
        index=["sample", "replicate"],
        columns="molecule",
        values="duration_sec"
    ).reset_index()

    cont_pivot.rename(columns={
        "dna": "containment_dna_sec",
        "protein": "containment_protein_sec"
    }, inplace=True)

    # --- Select relevant sketch timing columns ---
    sig_clean = sig_df[["sample", "replicate", "dna_runtime_sec", "protein_runtime_sec"]]

    # --- Merge all three ---
    merged = sig_clean.merge(cont_pivot, on=["sample", "replicate"], how="left")
    merged = merged.merge(fmh_df, on=["sample", "replicate"], how="left")

    # --- Compute total time ---
    merged["total_time_sec"] = (
        merged["dna_runtime_sec"].fillna(0)
        + merged["protein_runtime_sec"].fillna(0)
        + merged["containment_dna_sec"].fillna(0)
        + merged["containment_protein_sec"].fillna(0)
        + merged["duration_sec"].fillna(0)  # FMH time
    )

    # Automatically derive number of pairwise comparisons (sample size squared)
    merged["pairwise_count"] = merged["sample"].astype(int) ** 2

    # Normalized timing
    merged["time_per_pairwise"] = merged["total_time_sec"] / merged["pairwise_count"]


    # --- Save ---
    merged.to_csv(out_path, index=False)
    print(f"✅ Merged timing file written to: {out_path}")


if __name__ == "__main__":
    main()
