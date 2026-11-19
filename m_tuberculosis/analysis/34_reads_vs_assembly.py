"""Do read-derived and assembly-derived isolates give the same values? With the tool.

    usage:  python analysis/34_reads_vs_assembly.py [cores]

The read pipeline (steps 07-08) also processed 299 runs that do have an assembly. For those isolates
the whole coding genome is written twice, once with the substitutions found from the assembly
(work/fasta/all_cds/<sample>) and once with those found from the reads (here, named <run>), and the
tool (reference mode, scale factor 1) compares all of them at once. The same pair of isolates then
has one value from assemblies and one from reads.

out: results/reads_vs_assembly_tool/isolates.tsv   per isolate: substitutions found each way, shared
     results/reads_vs_assembly_tool/pairs.tsv      per pair of isolates: dN/dS both ways
     results/reads_vs_assembly_tool/summary.txt
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import anchored

MTB = anchored.MTB
TOOL = f"{MTB}/external/dnds-using-fmh/fmh_dnds.py"
WORK = f"{MTB}/work/reads_vs_assembly_tool"
OUT = f"{MTB}/results/reads_vs_assembly_tool"
CORES = sys.argv[1] if len(sys.argv) > 1 else "40"
os.makedirs(f"{WORK}/fasta", exist_ok=True)
os.makedirs(OUT, exist_ok=True)

calibration = pd.read_csv(f"{MTB}/data/isolates/calibration_runs.tsv", sep="\t")
both = [(run, sample) for run, sample in zip(calibration.run, calibration["sample"])
        if os.path.exists(f"{MTB}/work/snv_reads/{run}.calls.tsv.gz")
        and os.path.exists(f"{MTB}/work/fasta/all_cds/{sample}.ffn")]
print(f"{len(both)} isolates with both kinds of input", flush=True)

backbone, genes = anchored.read_backbone(), anchored.gene_table()
rows = []
for run, sample in both:
    from_reads = anchored.substitutions_of(run, "reads")
    from_assembly = anchored.substitutions_of(sample, "assembly")
    shared = set(from_reads.items()) & set(from_assembly.items())
    rows.append(dict(run=run, sample=sample, from_reads=len(from_reads), from_assembly=len(from_assembly),
                     shared=len(shared),
                     jaccard=len(shared) / max(len(set(from_reads.items()) | set(from_assembly.items())), 1)))
    with open(f"{WORK}/fasta/{run}.ffn", "w") as dna, open(f"{WORK}/fasta/{run}.faa", "w") as protein:
        for locus_tag, group, in_all, cds, translated in anchored.coding_sequences(anchored.apply(backbone, from_reads), genes):
            if in_all:
                dna.write(f">{locus_tag}\n{cds}\n")
                if translated:
                    protein.write(f">{locus_tag}\n{translated}\n")
isolates = pd.DataFrame(rows)
isolates.to_csv(f"{OUT}/isolates.tsv", sep="\t", index=False)

sheet = ["name,genome_filename,protein_filename",
         f"ancestral_reference,{MTB}/work/fasta/all_cds/ancestral_reference.ffn,{MTB}/work/fasta/all_cds/ancestral_reference.faa"]
sheet += [f"{sample},{MTB}/work/fasta/all_cds/{sample}.ffn,{MTB}/work/fasta/all_cds/{sample}.faa" for _, sample in both]
sheet += [f"{run},{WORK}/fasta/{run}.ffn,{WORK}/fasta/{run}.faa" for run, _ in both]
open(f"{WORK}/genomes.csv", "w").write("\n".join(sheet) + "\n")
environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
log = subprocess.run([sys.executable, TOOL, "--mode", "reference", "--reference", "ancestral_reference",
                      "--fasta_input_list", f"{WORK}/genomes.csv", "--scaled_input", "1", "--ksize", "7",
                      "--threshold", "0", "--cores", CORES, "--sketch_shards", "20", "--directory", f"{WORK}/tool"],
                     capture_output=True, text=True, env=environment, cwd=os.path.dirname(TOOL))
open(f"{WORK}/tool.log", "w").write(log.stdout + log.stderr)
if log.returncode != 0:
    sys.exit(f"the tool failed; see {WORK}/tool.log")

table = pd.read_parquet(f"{WORK}/tool/reference_scaled1/pairs")
value = {}
for a, b, d in zip(table.query_name.astype(str), table.match_name.astype(str), table["dN/dS"]):
    value[(a, b)] = value[(b, a)] = d
pairs = []
for index, (run_a, sample_a) in enumerate(both):
    for run_b, sample_b in both[index + 1:]:
        pairs.append(dict(sample_a=sample_a, sample_b=sample_b,
                          assembly_dnds=value[(sample_a, sample_b)], reads_dnds=value[(run_a, run_b)]))
pairs = pd.DataFrame(pairs).replace([np.inf, -np.inf], np.nan)
pairs.to_csv(f"{OUT}/pairs.tsv", sep="\t", index=False)
usable = pairs.dropna()
lines = [f"isolates with both kinds of input: {len(isolates)}",
         "substitutions per isolate, found from reads, from the assembly, and in both:",
         isolates[["from_reads", "from_assembly", "shared", "jaccard"]].describe().round(4).to_string(), "",
         f"pairs of isolates: {len(pairs)} ({len(usable)} with a finite value both ways)",
         f"median dN/dS from assemblies: {usable.assembly_dnds.median():.4f}",
         f"median dN/dS from reads: {usable.reads_dnds.median():.4f}",
         f"median difference, reads minus assemblies: {(usable.reads_dnds - usable.assembly_dnds).median():.4f}",
         f"Spearman correlation: {usable.assembly_dnds.corr(usable.reads_dnds, method='spearman'):.4f}",
         f"Pearson correlation: {usable.assembly_dnds.corr(usable.reads_dnds):.4f}"]
open(f"{OUT}/summary.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
