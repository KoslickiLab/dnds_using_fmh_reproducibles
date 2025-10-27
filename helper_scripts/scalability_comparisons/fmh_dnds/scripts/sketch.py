#!/usr/bin/env python3
import os
import subprocess
import argparse
import pandas as pd
import sys
import time
from datetime import datetime


def run_command(cmd):
    """Run a shell command safely with logging."""
    print(f"→ Running: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running command: {cmd}")
        sys.exit(e.returncode)

def main():
    parser = argparse.ArgumentParser(description="Sketch DNA and protein FASTA sequences using sourmash.")
    parser.add_argument("--dataset", required=True, help="CSV file with columns: name,dna_fasta,protein_fasta")
    parser.add_argument("--outdir", required=True, help="Base directory for sketches")
    parser.add_argument("--ksize", type=int, default=7, help="K-mer size (default: 7)")
    parser.add_argument("--scaled", type=int, default=1, help="Scaling factor (default: 1)")
    args = parser.parse_args()

    df = pd.read_csv(args.dataset)
    print(f"📂 Loaded {len(df)} entries from {args.dataset}")

    # Prepare timing log
    timing_log = []

    for idx, row in df.iterrows():
        name, dna_fasta, protein_fasta = row[0], row[1], row[2]

        sample_dir = os.path.dirname(dna_fasta)
        sig_dir = os.path.join(sample_dir, "signatures")
        os.makedirs(sig_dir, exist_ok=True)

        dna_sigfile = os.path.join(sig_dir, f"{name}_dna.sig")
        protein_sigfile = os.path.join(sig_dir, f"{name}_protein.sig")

        dna_cmd = f"sourmash sketch dna -p k={args.ksize},scaled={args.scaled} -o {dna_sigfile} {dna_fasta} --singleton"
        protein_cmd = f"sourmash sketch protein -p k={args.ksize},scaled={args.scaled} -o {protein_sigfile} {protein_fasta} --singleton"

        start_time = time.time()
        run_command(dna_cmd)
        dna_elapsed = time.time() - start_time


        start_time = time.time()
        run_command(protein_cmd)
        protein_elapsed = time.time() - start_time


        # Save timing info
        timing_log.append({
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "sample": name,
            "dna_fasta": dna_fasta,
            "protein_fasta": protein_fasta,
            "dna_runtime_sec": round(dna_elapsed, 2),
            "protein_runtime_sec": round(protein_elapsed, 2)
        })

    # Write timing results to CSV
    timing_df = pd.DataFrame(timing_log)
    timing_csv = os.path.join(args.outdir, "signature_generation_times.csv")
    timing_df.to_csv(timing_csv, index=False)


    print("✅ Finished generating all signatures.")
    print(f"🕒 Timing summary saved to: {timing_csv}")

if __name__ == "__main__":
    main()
