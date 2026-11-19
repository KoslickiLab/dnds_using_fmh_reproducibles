"""The input files of the FracMinHash dN/dS tool: one coding-sequence FASTA and one protein FASTA per
isolate and per group of genes, and the csv the tool reads.

    usage:  python analysis/30_write_inputs.py [--jobs N]

Isolates: the 39,379 with an AllTheBacteria assembly (substitutions from step 14) and the 10,824
without one whose reads were processed (steps 07-08). The whole coding genome (all_cds) is written
for all of them; the four gene groups for the assembly-derived isolates, which the group comparison
uses. The reference itself, with no substitutions, is written once per group under the name
`ancestral_reference`; it is the genome the tool's reference mode compares every isolate with.

out: work/fasta/<group>/<name>.ffn, <name>.faa      one record per coding sequence (locus tag)
     work/fasta/<group>/genomes.csv                 name,genome_filename,protein_filename
     work/fasta/assembly_names.txt                  the 39,379 assembly-derived isolates
     results/inputs/isolates.tsv                    per isolate: source and number of substitutions
"""
import argparse
import glob
import os
import sys
from multiprocessing import Pool

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import anchored

MTB = anchored.MTB
OUT = f"{MTB}/work/fasta"
REFERENCE_NAME = "ancestral_reference"


def isolates():
    """[(name, source)]: assembly-derived isolates, then read-derived ones without an assembly."""
    assembly = sorted(os.path.basename(p)[:-4] for p in glob.glob(f"{MTB}/work/variants/*.cov"))
    matched = pd.read_csv(f"{MTB}/data/isolates/gan_to_atb.tsv", sep="\t", dtype=str)
    run_to_sample = dict(zip(matched.run, matched["sample"].fillna("")))
    have = set(assembly)
    reads = sorted(os.path.basename(p)[:-len(".calls.tsv.gz")]
                   for p in glob.glob(f"{MTB}/work/snv_reads/*.calls.tsv.gz"))
    reads = [run for run in reads if run_to_sample.get(run, "") not in have]
    return [(name, "assembly") for name in assembly] + [(name, "reads") for name in reads]


def write(task):
    name, source = task
    groups = anchored.GROUPS if source in ("assembly", "reference") else ["all_cds"]
    if all(os.path.exists(f"{OUT}/{g}/{name}.faa") for g in groups):
        return name, source, None
    substitutions = {} if source == "reference" else anchored.substitutions_of(name, source)
    sequences = anchored.coding_sequences(anchored.apply(BACKBONE, substitutions), GENES)
    for group in groups:
        with open(f"{OUT}/{group}/{name}.ffn.tmp", "w") as dna, open(f"{OUT}/{group}/{name}.faa.tmp", "w") as protein:
            for locus_tag, gene_group, in_all, cds, translated in sequences:
                if anchored.members(gene_group, in_all, group):
                    dna.write(f">{locus_tag}\n{cds}\n")
                    if translated:
                        protein.write(f">{locus_tag}\n{translated}\n")
        os.replace(f"{OUT}/{group}/{name}.ffn.tmp", f"{OUT}/{group}/{name}.ffn")
        os.replace(f"{OUT}/{group}/{name}.faa.tmp", f"{OUT}/{group}/{name}.faa")
    return name, source, len(substitutions)


def initialise():
    global BACKBONE, GENES
    BACKBONE = anchored.read_backbone()
    GENES = anchored.gene_table()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs", type=int, default=200)
    arguments = parser.parse_args()
    for group in anchored.GROUPS:
        os.makedirs(f"{OUT}/{group}", exist_ok=True)
    os.makedirs(f"{MTB}/results/inputs", exist_ok=True)

    tasks = [(REFERENCE_NAME, "reference")] + isolates()
    print(f"{len(tasks) - 1} isolates", flush=True)
    rows = []
    with Pool(arguments.jobs, initializer=initialise) as pool:
        for i, row in enumerate(pool.imap_unordered(write, tasks, chunksize=8), 1):
            rows.append(row)
            if i % 5000 == 0:
                print(f"  {i}/{len(tasks)}", flush=True)

    table = pd.DataFrame(rows, columns=["name", "source", "substitutions"])
    previous = f"{MTB}/results/inputs/isolates.tsv"
    if os.path.exists(previous):            # isolates written by an earlier run keep their counts
        earlier = pd.read_csv(previous, sep="\t").set_index("name").substitutions
        table["substitutions"] = table["substitutions"].fillna(table["name"].map(earlier))
    table.sort_values(["source", "name"]).to_csv(previous, sep="\t", index=False)

    for group in anchored.GROUPS:
        names = [name for name, source, _ in rows if group == "all_cds" or source in ("assembly", "reference")]
        lines = ["name,genome_filename,protein_filename"] + [
            f"{name},{OUT}/{group}/{name}.ffn,{OUT}/{group}/{name}.faa" for name in sorted(names)]
        open(f"{OUT}/{group}/genomes.csv", "w").write("\n".join(lines) + "\n")
        print(f"{group}: {len(names)} genomes (including the reference)")
    assembly = sorted(name for name, source, _ in rows if source == "assembly")
    open(f"{OUT}/assembly_names.txt", "w").write("\n".join(assembly) + "\n")
    print(table.groupby("source").substitutions.describe().to_string())
