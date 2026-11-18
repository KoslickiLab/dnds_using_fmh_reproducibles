#!/usr/bin/env python
"""Check that the analysis starts from the same data as the manuscript.

    usage:  python checks/check_against_gtdb.py

1. The Prodigal gene calls made in step 2 are compared with the gene calls of GTDB release 214
   for the 243 genomes of the hierarchical edge bundling figure
   (checks/gtdb_r214_gene_calls_md5.tsv holds the md5 checksums of the GTDB protein and CDS
   sequences of each genome, taken from gtdb_proteins_aa_reps_r214 and gtdb_proteins_nt_reps_r214).
2. The FracMinHash dN/dS computed in step 4 is compared with the value of each pair in the table
   behind the manuscript's figure (all_pairs/pairs.tsv), for the 140 pairs (step 11) and for the
   pairs examined in detail (pairs/examined_pairs.tsv).

Known difference: in GCA_021812785.1, the start codon of one of 2,488 genes is called 9 bp
further downstream than in GTDB; the FracMinHash dN/dS of its pair is unaffected.
"""
import hashlib
import os
import pandas as pd


def md5_of_sequences(path, strip_star):
    """md5 of all sequence lines of a FASTA file, joined, upper case (stop symbols '*' removed for proteins)."""
    digest = hashlib.md5()
    for line in open(path):
        if not line.startswith(">"):
            sequence = line.strip()
            if strip_star:
                sequence = sequence.replace("*", "")
            digest.update(sequence.upper().encode())
    return digest.hexdigest()


# ---- 1. gene calls
different, checked = [], 0
for row in pd.read_csv("checks/gtdb_r214_gene_calls_md5.tsv", sep="\t").itertuples():
    if not os.path.exists(f"data/genes/{row.accession}.faa"):
        continue                    # genome not downloaded (the 140-pair analysis was not run)
    checked += 1
    same = (md5_of_sequences(f"data/genes/{row.accession}.faa", True) == row.proteins_md5
            and md5_of_sequences(f"data/genes/{row.accession}.ffn", False) == row.cds_md5)
    if not same:
        different.append(row.accession)
print(f"gene calls: {checked - len(different)} of {checked} genomes identical to GTDB R214"
      + (f"; different: {', '.join(different)}" if different else ""))

# ---- 2. FracMinHash dN/dS
paper = pd.read_csv("all_pairs/pairs.tsv", sep="\t").set_index("pair").fmh_dnds_in_paper
folders = {f"work/all_pairs/{pair}": pair for pair in paper.index}
examined = pd.read_csv("pairs/examined_pairs.tsv", sep="\t")
folders.update(dict(zip(examined.folder, examined.pair)))
largest, compared = 0.0, 0
for folder, pair in folders.items():
    if not os.path.exists(f"{folder}/results/fmh_dnds.tsv"):
        continue
    ours = pd.read_csv(f"{folder}/results/fmh_dnds.tsv", sep="\t").set_index("genes").loc["all", "fmh_dnds"]
    compared += 1
    largest = max(largest, abs(ours - paper[pair]))
    if abs(ours - paper[pair]) > 1e-6:
        print(f"{folder}: FracMinHash dN/dS {ours:.6f}, manuscript {paper[pair]:.6f} (DIFFERENT)")
print(f"FracMinHash dN/dS: {compared} comparisons, largest difference to the manuscript's value {largest:.1e}")

all_good = different in ([], ["GCA_021812785.1"]) and largest < 1e-6
print("all checks passed" if all_good else "SOME CHECKS FAILED")
