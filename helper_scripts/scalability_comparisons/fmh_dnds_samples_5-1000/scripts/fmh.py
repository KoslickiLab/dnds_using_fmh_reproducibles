#!/usr/bin/env python3
import os
import subprocess
import time
import pandas as pd
import argparse


def run_command(cmd):
    """Run a command and return elapsed time in seconds."""
    print(f"→ Running: {cmd}")
    start = time.time()
    try:
        subprocess.run(cmd, shell=True, check=True)
        status = "success"
    except subprocess.CalledProcessError as e:
        print(f"❌ Command failed: {e}")
        status = "failed"
    duration = round(time.time() - start, 5)
    print(f"✓ Completed in {duration} seconds ({status})")
    return duration, status


def main():
    parser = argparse.ArgumentParser(description="Run dN/dS scalability comparison workflow.")
    parser.add_argument("--base",required=True,help="Base directory for the scalability comparison (e.g. /path/to/base)")
    parser.add_argument("--samples",nargs="+",type=int,required=True,help="List of sample sizes to process (e.g. 5 10 100 1000)")
    parser.add_argument("--results_dir",required=True,help="Directory where containments can be found")
    parser.add_argument("--dataset_csv",required=True,help="Path to dataset.csv file")
    parser.add_argument("--ksize",type=int,default=15,help="K-mer size (default: 15)")
    parser.add_argument("--out_log",required=True,help="Output CSV file for logging results")
    args = parser.parse_args()


    df = pd.read_csv(args.dataset_csv)
    results = []

    # dataframe contains column: name = sample5_cds_WP_000551270.1
    for idx, row in df.iterrows():
        full_name = row["name"]

        # example full_name → "sample5_cds_WP_000551270.1"
        try:
            wp_id = full_name.split("_cds_")[1]
        except:
            print(f"⚠️ Could not parse WP id from: {full_name}")
            continue

        # deduce sample number
        sample_part = full_name.split("_cds_")[0]   # "sample5"
        s = sample_part.replace("sample", "")       # "5"

        cont_dir = os.path.join(args.base, args.results_dir, f"sample_{s}", "containments")

        dna_cmp = os.path.join(cont_dir, f"compare_{wp_id}.dna.{args.ksize}.csv")
        prot_cmp = os.path.join(cont_dir, f"compare_{wp_id}.protein.{args.ksize}.csv")

        fmh_dir = os.path.join(args.base, args.results_dir, f"sample_{s}")

        if not os.path.exists(dna_cmp):
            print(f"❌ Missing DNA compare file: {dna_cmp}")
            continue
        if not os.path.exists(prot_cmp):
            print(f"❌ Missing protein compare file: {prot_cmp}")
            continue

        cmd = (
            f"python3 /data/jzr5814/repositories/dnds-using-fmh/src/script_fmh_omega.py "
            f"--dna_sample {dna_cmp} --protein_sample {prot_cmp} "
            f"--scaled_input 1 --ksize {args.ksize} --mode test "
            f"--directory {fmh_dir} --cores 10 --threshold 0"
        )

        print(f"\n=== Running {full_name}")
        duration, status = run_command(cmd)

        results.append({
            "sample": s,
            "wp_id": wp_id,
            "duration_sec": duration,
            "status": status
        })

    pd.DataFrame(results).to_csv(args.out_log, index=False)
    print(f"\n🕒 Timing results saved to {args.out_log}")


if __name__ == "__main__":
    main()
