#!/usr/bin/env python3
import argparse
from pathlib import Path
import subprocess
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", nargs="+", required=True, type=int, help="Sample IDs")
    parser.add_argument("--replicates", type=int, required=True, help="Number of replicates per sample")
    parser.add_argument("--base", required=True, help="Base directory for all samples")
    parser.add_argument("--translator", required=True, help="Path to translation program (transeq)")
    parser.add_argument("--dataset", required=True, help="Output dataset CSV path")
    args = parser.parse_args()

    base = Path(args.base)
    rows = []

    for s in args.samples:
        for r in range(1, args.replicates + 1):
            fna_file = base / f"sample_{s}" / f"rep_{r}" / "sampled_cds.fna"
            protein_file = base / f"sample_{s}" / f"rep_{r}" / "translated_protein.faa"

            # Make sure directory exists
            protein_file.parent.mkdir(parents=True, exist_ok=True)

            # Translate if protein file does not exist
            if not protein_file.exists():
                subprocess.run([args.translator, "-sequence", str(fna_file), "-outseq", str(protein_file)],
                               check=True)

            # Add a row to dataset
            rows.append({
                "name": f"sample{s}_rep{r}",
                "genome_filename": str(fna_file.resolve()),
                "protein_filename": str(protein_file.resolve())
            })

    # Save dataset CSV
    df = pd.DataFrame(rows)
    df.to_csv(args.dataset, index=False)
    print(f"Dataset CSV written to {args.dataset}")

if __name__ == "__main__":
    main()
