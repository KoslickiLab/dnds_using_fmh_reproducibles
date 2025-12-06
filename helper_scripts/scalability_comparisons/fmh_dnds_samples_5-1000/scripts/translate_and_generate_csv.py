#!/usr/bin/env python3
import argparse
from pathlib import Path
import subprocess
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", nargs="+", required=True, type=int, help="Sample IDs")
    parser.add_argument("--base", required=True, help="Base directory for all samples")
    parser.add_argument("--translator", required=True, help="Path to translation program (transeq)")
    parser.add_argument("--dataset", required=True, help="Output dataset CSV path")
    args = parser.parse_args()

    base = Path(args.base)
    rows = []

    for s in args.samples:
        sample_dir = base /  f"results_from_script_testing/sample_{s}" / f"fasta"   # input CDS location

        # NEW: output location for translated files
        out_dir = base / "snakemake_fmh_comparison_samples_10-1000" / f"sample_{s}" / "faa"
        out_dir.mkdir(parents=True, exist_ok=True)

        # Find ALL CDS files
        fna_files = list(sample_dir.glob("cds_*.fasta"))
        if not fna_files:
            raise FileNotFoundError(f"No CDS fasta/fna files found in {sample_dir}")

        for fna_file in fna_files:

            # Build protein output filename inside the new directory
            protein_file = out_dir / (fna_file.stem + "_translated_protein.faa")

            # Translate if needed
            if not protein_file.exists():
                subprocess.run(
                    [args.translator, "-sequence", str(fna_file), "-outseq", str(protein_file)],
                    check=True
                )

            # Add row PER FILE
            rows.append({
                "name": f"sample{s}_{fna_file.stem}",
                "genome_filename": str(fna_file.resolve()),
                "protein_filename": str(protein_file.resolve())
            })


    # Save dataset CSV
    df = pd.DataFrame(rows)
    df.to_csv(args.dataset, index=False)
    print(f"Dataset CSV written to {args.dataset}")

if __name__ == "__main__":
    main()
