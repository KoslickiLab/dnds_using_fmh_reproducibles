#!/usr/bin/env python
"""Step 14. The one pair of related genomes with a FracMinHash dN/dS above 1, and the
translation table of its gene calls.

    usage:  python analysis/14_translation_table_check.py      (after analysis/11_all_pairs_run.sh)

Among the 140 pairs, one pair with an estimate above 1 has no near-identical DNA to explain it:
Neisseria lactamica GCA_000193795.2 and Neisseria meningitidis GCA_923184405.1 (pair P128).
GCA_923184405.1 has no NCBI annotation, and by GTDB's rule (step 2) its genes were predicted
with translation table 4, in which TGA codes for tryptophan. Neisseria uses table 11. Here the
genes of that genome are predicted again with table 11 and the FracMinHash dN/dS of the pair is
recomputed, to show how much of the estimate is due to the gene calls.

Output: all_pairs/results/translation_table_check.tsv
        (translation_table, genes, median_protein_length, fmh_dnds of the pair)
"""
import os
import subprocess
import sys
import pandas as pd
from common import read_fasta, THREADS

PARTNER, GENOME = "GCA_000193795.2", "GCA_923184405.1"
work = os.path.abspath("work/translation_table_check")
os.makedirs(work, exist_ok=True)

# gene calls of GCA_923184405.1 with table 11 (step 2 kept the table-4 calls in data/genes/)
subprocess.run(["prodigal", "-m", "-p", "single", "-g", "11", "-q", "-i", f"data/genomes/{GENOME}/genome.fna",
                "-d", f"{work}/table11.ffn", "-a", f"{work}/table11.faa", "-o", os.devnull], check=True)
genes = {4: f"data/genes/{GENOME}", 11: f"{work}/table11"}

sample_sheet = ["name,genome_filename,protein_filename",
                f"partner,{os.path.abspath(f'data/genes/{PARTNER}.ffn')},{os.path.abspath(f'data/genes/{PARTNER}.faa')}"]
for table, prefix in genes.items():
    sample_sheet.append(f"table{table},{os.path.abspath(prefix + '.ffn')},{os.path.abspath(prefix + '.faa')}")
open(f"{work}/samples.csv", "w").write("\n".join(sample_sheet) + "\n")

# FracMinHash dN/dS with the parameters of step 4
for old in ("fmh_omega_7.csv", "data.zip"):
    if os.path.exists(f"{work}/{old}"):
        os.remove(f"{work}/{old}")
command = ["python", "external/dnds-using-fmh/fmh_dnds.py", "--fasta_input_list", f"{work}/samples.csv",
           "--ksize", "7", "--scaled_input", "500", "--directory", work, "--cores", THREADS, "--threshold", "0"]
log = subprocess.run(command, capture_output=True, text=True)
open(f"{work}/fmh_dnds.log", "w").write(log.stdout + log.stderr)
if log.returncode != 0:
    sys.exit(f"fmh_dnds.py failed; see {work}/fmh_dnds.log")
results = pd.read_csv(f"{work}/fmh_omega_7.csv")

rows = []
for table, prefix in genes.items():
    proteins = read_fasta(prefix + ".faa")
    lengths = pd.Series([len(protein.rstrip("*")) for protein in proteins.values()])
    wanted = {"partner", f"table{table}"}
    value = results[[{q, m} == wanted for q, m in zip(results.query_name, results.match_name)]].iloc[0]["dN/dS"]
    rows.append(dict(translation_table=table, genes=len(proteins), median_protein_length=lengths.median(), fmh_dnds=value))
table = pd.DataFrame(rows)
table.to_csv("all_pairs/results/translation_table_check.tsv", sep="\t", index=False)
print(table.to_string(index=False))
