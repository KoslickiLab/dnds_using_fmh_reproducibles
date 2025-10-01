#!/usr/bin/env python3

import argparse
import random
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Randomly sample X genes from a FASTA file.")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file with sampled genes")
    parser.add_argument("--sample_n", type=int, required=True, help="Number of genes to sample")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")

    args = parser.parse_args()

    # Set seed for reproducibility
    if args.seed is not None:
        random.seed(args.seed)

    # Load sequences
    records = list(SeqIO.parse(args.input, "fasta"))
    if len(records) == 0:
        raise ValueError(f"No records found in {args.input}")

    if args.sample_n > len(records):
        raise ValueError(
            f"Requested {args.sample_n} sequences, but only {len(records)} available."
        )

    # Randomly sample
    sampled_records = random.sample(records, args.sample_n)

    # Write output
    SeqIO.write(sampled_records, args.output, "fasta")

    print(
        f"Sampled {args.sample_n} genes from {args.input} "
        f"and wrote them to {args.output}"
    )

if __name__ == "__main__":
    main()
