#!/usr/bin/env python3
import os
import argparse
import pandas as pd

def extract_sample_number(df, sample_col):
    """Extract numeric sample from the sample column."""
    df["sample"] = df[sample_col].str.extract(r"sample(\d+)", expand=False).astype(int)
    return df

def main():
    parser = argparse.ArgumentParser(description="Merge all timing logs into one dataset using wp_id.")
    parser.add_argument("--base", required=True, help="Base directory containing timing logs.")
    args = parser.parse_args()

    base = args.base

    # --- File paths ---
    sig_path = os.path.join(base, "signature_generation_times.csv")
    cont_path = os.path.join(base, "containment_times.csv")
    fmh_path = os.path.join(base, "fmh_omega_times.csv")
    out_path = os.path.join(base, "all_timings_merged.csv")

    # --- Load CSVs ---
    sig_df = pd.read_csv(sig_path)
    cont_df = pd.read_csv(cont_path)
    fmh_df = pd.read_csv(fmh_path)

    # --- Extract numeric sample and wp_id ---
    sig_df["wp_id"] = sig_df["sample"].apply(lambda x: "_".join(x.split("_")[-2:]))
    sig_df["sample"] = sig_df["sample"].str.extract(r"sample(\d+)", expand=False).astype(int)
    
    # --- Pivot containment so DNA & protein times become columns ---
    cont_pivot = cont_df.pivot_table(
        index="wp_id",
        columns="molecule",
        values="duration_sec"
    ).reset_index()
    cont_pivot.rename(columns={
        "dna": "containment_dna_sec",
        "protein": "containment_protein_sec"
    }, inplace=True)

    # --- Select relevant signature timing columns ---
    print(sig_df.columns)
    sig_clean = sig_df[["sample", "wp_id", "dna_runtime_sec", "protein_runtime_sec"]]

    # --- Merge signature times with containment and FMH times using wp_id ---
    merged = sig_clean.merge(cont_pivot, on="wp_id", how="left")
    merged = merged.merge(fmh_df[["wp_id", "duration_sec"]], on="wp_id", how="left")

    # --- Compute total time ---
    merged["total_time_sec"] = (
        merged["dna_runtime_sec"].fillna(0)
        + merged["protein_runtime_sec"].fillna(0)
        + merged.get("containment_dna_sec", 0).fillna(0)
        + merged.get("containment_protein_sec", 0).fillna(0)
        + merged.get("duration_sec", 0).fillna(0)  # FMH time
    )

    # --- Save merged CSV ---
    merged.to_csv(out_path, index=False)
    print(f"✅ Merged timing file written to: {out_path}")

if __name__ == "__main__":
    main()
