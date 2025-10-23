#!/usr/bin/env python3
import os
import subprocess
import argparse
import pandas as pd
import sys

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

    for idx, row in df.iterrows():
        name, dna_fasta, protein_fasta = row[0], row[1], row[2]

        sample_dir = os.path.dirname(dna_fasta)
        sig_dir = os.path.join(sample_dir, "signatures")
        os.makedirs(sig_dir, exist_ok=True)

        dna_sigfile = os.path.join(sig_dir, f"{name}_dna.sig")
        protein_sigfile = os.path.join(sig_dir, f"{name}_protein.sig")

        dna_cmd = f"sourmash sketch dna -p k={args.ksize},scaled={args.scaled} -o {dna_sigfile} {dna_fasta} --singleton"
        protein_cmd = f"sourmash sketch protein -p k={args.ksize},scaled={args.scaled} -o {protein_sigfile} {protein_fasta} --singleton"

        run_command(dna_cmd)
        run_command(protein_cmd)

        # --- Protein sketch ---
        protein_sigfile = os.path.join(sig_dir, f"{name}_protein.sig")
        protein_cmd = (
            f"sourmash sketch protein -p k={args.ksize},scaled={args.scaled} "
            f"-o {protein_sigfile} {protein_fasta} --singleton"
        )
        run_command(protein_cmd)

    print("✅ Finished generating all signatures.")

if __name__ == "__main__":
    main()
