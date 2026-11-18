#!/usr/bin/env python
"""Step 13. Experiment: does an incomplete overlap of gene content raise FracMinHash dN/dS?

    usage:  python analysis/13_gene_removal_experiment.py      (after the pairs faecalibacterium_roseburia and
                                                                roseburia have been run, which download the genomes)

Two complete genomes of ONE species diverge little, so their FracMinHash dN/dS is well below 1.
Here we delete a random set of genes from each of the two genomes (independently), which mimics
genomes that are incomplete (MAGs) or that differ in gene content, and recompute FracMinHash
dN/dS. Nothing about the divergence of the remaining genes changes.

Two pairs of conspecific genomes that are among the genomes of the pairs examined in detail are used.
Each genome keeps 100, 90, 80, 70, 60 or 50% of its genes; every setting is repeated five times
with different random genes.

Expectation: if each genome keeps a fraction p of its genes, a k-mer of one genome still has its
counterpart in the other with probability p, so both containments shrink by the factor p:
        C_nt(p) = p * C_nt(1),   C_aa(p) = p * C_aa(1).

Output: all_pairs/results/gene_removal.tsv
        (pair, kept_percent, replicate, containment_nt, containment_aa, fmh_dnds)
"""
import os
import random
import subprocess
import sys

import pandas as pd
from common import read_fasta, THREADS

PAIRS = {    # name -> the two genome accessions (both complete or chromosome-level, same species)
    "F. duncaniae A2-165 vs NB2B-19-DCM": ("GCF_010509575.1", "GCF_024167885.1"),
    "R. intestinalis L1-82 vs NB2A 31_NB": ("GCF_900537995.1", "GCF_023497585.1"),
}
KEPT_PERCENT = [100, 90, 80, 70, 60, 50]
REPLICATES = 5

out_dir = "all_pairs/results"
work = os.path.abspath("work/gene_removal")
os.makedirs(work, exist_ok=True)

# 1. write the reduced genomes
random.seed(2026)
sample_sheet = ["name,genome_filename,protein_filename"]
runs = []      # (pair name, kept percent, replicate, label of genome 1, label of genome 2)
for number, (pair, accessions) in enumerate(PAIRS.items(), start=1):
    for kept in KEPT_PERCENT:
        for replicate in range(1, REPLICATES + 1 if kept < 100 else 2):       # 100% needs no repeats
            labels = []
            for side, accession in zip("AB", accessions):
                cds = read_fasta(f"data/genes/{accession}.ffn")
                proteins = read_fasta(f"data/genes/{accession}.faa")
                genes = list(cds)
                keep = set(random.sample(genes, round(len(genes) * kept / 100)))
                label = f"pair{number}_{side}_keep{kept}_rep{replicate}"
                with open(f"{work}/{label}.ffn", "w") as out:
                    out.write("".join(f">{g}\n{cds[g]}\n" for g in genes if g in keep))
                with open(f"{work}/{label}.faa", "w") as out:
                    out.write("".join(f">{g}\n{proteins[g]}\n" for g in genes if g in keep))
                sample_sheet.append(f"{label},{work}/{label}.ffn,{work}/{label}.faa")
                labels.append(label)
            runs.append((pair, kept, replicate, labels[0], labels[1]))
open(f"{work}/samples.csv", "w").write("\n".join(sample_sheet) + "\n")

# 2. FracMinHash dN/dS (same parameters as in step 4)
for old in ("fmh_omega_7.csv", "data.zip"):
    if os.path.exists(f"{work}/{old}"):
        os.remove(f"{work}/{old}")
command = ["python", "external/dnds-using-fmh/fmh_dnds.py", "--fasta_input_list", f"{work}/samples.csv",
           "--ksize", "7", "--scaled_input", "500", "--directory", work, "--cores", THREADS, "--threshold", "0"]
log = subprocess.run(command, capture_output=True, text=True)
open(f"{work}/fmh_dnds.log", "w").write(log.stdout + log.stderr)
if log.returncode != 0:
    sys.exit(f"fmh_dnds.py failed; see {work}/fmh_dnds.log")

# 3. collect the comparisons of interest
results = pd.read_csv(f"{work}/fmh_omega_7.csv")
by_names = {frozenset((q, m)): (c_nt, c_aa, value) for q, m, c_nt, c_aa, value in
            zip(results.query_name, results.match_name, results.DNA_max_Cfrac, results.AA_max_Cfrac, results["dN/dS"])}
rows = []
for pair, kept, replicate, label_a, label_b in runs:
    c_nt, c_aa, value = by_names[frozenset((label_a, label_b))]
    rows.append(dict(pair=pair, kept_percent=kept, replicate=replicate, containment_nt=c_nt, containment_aa=c_aa, fmh_dnds=value))
table = pd.DataFrame(rows)
table.to_csv(f"{out_dir}/gene_removal.tsv", sep="\t", index=False)
print(table.groupby(["pair", "kept_percent"]).fmh_dnds.agg(["mean", "min", "max"]).round(2).to_string())
