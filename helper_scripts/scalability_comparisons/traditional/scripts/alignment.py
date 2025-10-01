#!/usr/bin/env python3
import argparse
from pathlib import Path
from Bio import SeqIO
import subprocess
import time
import csv
import itertools

def write_pair_fasta(seq1, seq2, out_file):
    """Write a temporary FASTA file for a sequence pair, sequences on one line."""
    with open(out_file, "w") as f:
        f.write(f">{seq1.id}\n{str(seq1.seq)}\n")
        f.write(f">{seq2.id}\n{str(seq2.seq)}\n")

def run_pairwise_clustal(fasta_file, outdir):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    timing_file = outdir / "timing_alignments.tsv"
    with open(timing_file, "w", newline="") as tf:
        writer = csv.writer(tf, delimiter="\t")
        writer.writerow(["seq1", "seq2", "time_sec", "aln_file"])
        
        # Iterate over all sequence pairs
        for rec1, rec2 in itertools.combinations(records, 2):
            aln_name = f"{rec1.id}_{rec2.id}.aln"
            aln_path = outdir / aln_name
            
            # Temporary pairwise FASTA file
            pair_fasta = outdir / f"{rec1.id}_{rec2.id}.fasta"
            write_pair_fasta(rec1, rec2, pair_fasta)
            
            start = time.time()
            try:
                # Run ClustalW on this pair
                subprocess.run(
                    ["clustalw2", "-INFILE=" + str(pair_fasta), "-OUTFILE=" + str(aln_path), "-QUIET"],
                    check=True
                )
            except subprocess.CalledProcessError:
                print(f"ClustalW failed for {rec1.id} vs {rec2.id}")
                continue
            elapsed = time.time() - start
            
            # Record timing
            writer.writerow([rec1.id, rec2.id, f"{elapsed:.4f}", str(aln_path)])
            
            # Clean up temporary pairwise FASTA
            pair_fasta.unlink()
    
    print(f"Done: {len(records)} sequences, timings written to {timing_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pairwise ClustalW alignments for sampled FASTA")
    parser.add_argument("--input", required=True, help="Sampled FASTA file")
    parser.add_argument("--out_dir", required=True, help="Output directory for alignments and timing")
    args = parser.parse_args()
    
    run_pairwise_clustal(args.input, args.out_dir)
