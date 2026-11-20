"""What the equivalence check of step 31 found, as tables.

    usage:  python checks/31b_check_summary.py

The log of step 31 holds every run of the check; the last run of each group is the one on the final
code. Every comparison must report identical values in every column.

out: results/checks/reference_mode.tsv    per group and scale factor: pairs compared, identical or not
     results/checks/main_mode_sketch.tsv  size of main mode's single sketch file for the whole coding
                                          genomes of the check (scale factor 1), and its size per genome
"""
import os
import re

import pandas as pd

MTB = os.getcwd()
OUT = f"{MTB}/results/checks"
os.makedirs(OUT, exist_ok=True)

last_run = {}
for block in open(f"{MTB}/logs/31_check_reference_mode.log").read().split("\n--- ")[1:]:
    last_run[block.split(":", 1)[0]] = block
rows = []
for group, block in last_run.items():
    for scale, pairs, verdict in re.findall(r"scale factor (\d+)[^:]*: (\d+) pairs, ([^\n]*)", block):
        rows.append(dict(group=group, scale=int(scale), pairs=int(pairs),
                         identical=verdict == "identical in every column", report=verdict))
table = pd.DataFrame(rows)
table.to_csv(f"{OUT}/reference_mode.tsv", sep="\t", index=False)
assert table.identical.all(), "the check found differences"
print(table.to_string(index=False))

genomes = sum(1 for _ in open(f"{MTB}/work/equivalence/all_cds_80.csv")) - 1
size = os.path.getsize(f"{MTB}/work/equivalence/all_cds/main/data.zip")
pd.DataFrame([dict(genomes=genomes, bytes=size, bytes_per_genome=size / genomes)]).to_csv(
    f"{OUT}/main_mode_sketch.tsv", sep="\t", index=False)
print(f"main mode's sketch file: {size:,} bytes for {genomes} whole coding genomes")
