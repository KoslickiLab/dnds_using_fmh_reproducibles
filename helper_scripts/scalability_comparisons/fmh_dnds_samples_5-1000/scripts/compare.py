#!/usr/bin/env python3
import os
import subprocess
import argparse
import pandas as pd
import time
import glob

def run_command(cmd):
    """Run a shell command and return its execution time in seconds."""
    start_time = time.time()
    subprocess.run(cmd, shell=True, check=True)
    return time.time() - start_time

def main():
    parser = argparse.ArgumentParser(description="Compare DNA and protein signatures with sourmash.")
    parser.add_argument("--base", required=True, help="Base directory with sample_X folders.")
    parser.add_argument("--ksize", type=int, default=7, help="K-mer size (default: 7)")
    parser.add_argument("--outcsv", required=True, help="Output CSV for all timing records")
    args = parser.parse_args()

    samples = [5, 10, 100, 1000]  # sample numbers
    all_records = []

    for s in samples:
        base_dir = os.path.join(args.base, f"sample_{s}")
        sig_dir = os.path.join(base_dir, "signatures")
        if not os.path.exists(sig_dir):
            print(f"[WARNING] Signatures directory not found for sample{s}, skipping.")
            continue

        # Find all DNA signature files for this sample
        dna_files = sorted(glob.glob(os.path.join(sig_dir, f"sample{s}_cds_*_dna.sig")))
        if not dna_files:
            print(f"[WARNING] No DNA signature files found for sample{s}, skipping.")
            continue

        for dna_ref in dna_files:
            # Automatically find the corresponding protein file
            prot_ref = dna_ref.replace("_dna.sig", "_protein.sig")
            if not os.path.exists(prot_ref):
                print(f"[WARNING] Matching protein signature not found for {dna_ref}, skipping this pair.")
                continue

            # Extract WP_ID from filename
            # e.g., 'sample5_cds_WP_000040453.1_dna.sig' -> 'WP_000040453.1'
            wp_id = os.path.basename(dna_ref).split("_cds_")[1].split("_dna")[0]

            # Loop over DNA and protein
            for molecule, ref in [("dna", dna_ref), ("protein", prot_ref)]:
                out_csv = os.path.join(
                    base_dir,
                    "containments",
                    f"compare_{wp_id}.{molecule}.{args.ksize}.csv"
                )
                os.makedirs(os.path.dirname(out_csv), exist_ok=True)

                cmd = f"sourmash compare {ref} {ref} --containment --{molecule} --ksize {args.ksize} --csv {out_csv}"
                print(f"Running: {cmd}")
                duration = run_command(cmd)

                all_records.append({
                    "sample": s,
                    "wp_id": wp_id,
                    "molecule": molecule,
                    "ksize": args.ksize,
                    "duration_sec": round(duration, 2)
                })

    # Save all timing records to the consolidated CSV
    df = pd.DataFrame(all_records)
    df.to_csv(args.outcsv, index=False)
    print(f"✅ All comparison timings saved to {args.outcsv}")

if __name__ == "__main__":
    main()
