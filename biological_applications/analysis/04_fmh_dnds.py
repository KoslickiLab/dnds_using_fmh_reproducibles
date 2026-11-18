#!/usr/bin/env python
"""Step 4. FracMinHash dN/dS of a pair, with and without the near-identical genes.

    usage:  python analysis/04_fmh_dnds.py <pair folder>

The FracMinHash dN/dS implementation (external/dnds-using-fmh/fmh_dnds.py, unmodified) is run on
several versions of the two genomes: A and B as they are, and A and B after deleting the genes
that step 3 called near-identical (at least 50% of their 21-mers shared). To show that the
result does not depend on that 50%, the deletion is repeated with cut-offs of 25% and 75%.
Parameters are those of the manuscript: k_aa = 7 (so k_nt = 21) and a scale factor of 500.
The containment threshold is set to 0 so that the pair is reported whatever its containment.
Finally the complete genomes are compared once more with a scale factor of 1, that is, with
every k-mer instead of a sketch, to measure how much the sketch itself contributes.

Output:
    <pair folder>/results/fmh_dnds.tsv   one row per version, with the nucleotide and protein
        containment and the FracMinHash dN/dS. Versions (column "genes"):
            all                                  complete genomes, scale factor 500 (the manuscript's setting)
            without_near_identical               near-identical genes removed from both genomes
            without_genes_sharing_25_percent     the same with a 25% cut-off
            without_genes_sharing_75_percent     the same with a 75% cut-off
            all_every_kmer                       complete genomes, scale factor 1
"""
import sys
import os
import subprocess
import pandas as pd
from common import pair_accessions, read_fasta, write_fasta, work_dir, THREADS

pair_dir = sys.argv[1]
work = work_dir(pair_dir, "fmh")
accession_a, accession_b = pair_accessions(pair_dir)

# name of the row in the output -> genes to delete (None = delete nothing, otherwise the minimum
# fraction of shared 21-mers at which a gene is deleted)
VERSIONS = {"all": None, "without_near_identical": 0.50,
            "without_genes_sharing_25_percent": 0.25, "without_genes_sharing_75_percent": 0.75}
MIN_KMERS = 100      # as in step 3: genes with fewer 21-mers are never deleted

# 1. write the inputs and the sample sheet that fmh_dnds.py expects
sample_sheet = ["name,genome_filename,protein_filename"]
for side, accession in (("A", accession_a), ("B", accession_b)):
    cds = read_fasta(f"data/genes/{accession}.ffn")
    proteins = read_fasta(f"data/genes/{accession}.faa")
    genes = pd.read_csv(f"{pair_dir}/results/genes_{side}.tsv", sep="\t")
    for version, cutoff in VERSIONS.items():
        if cutoff is None:
            drop = set()
        else:
            drop = set(genes[(genes.fraction_shared >= cutoff) & (genes.n_kmers >= MIN_KMERS)].gene)
        label = f"{side}_{version}"
        write_fasta({g: s for g, s in cds.items() if g not in drop}, f"{work}/{label}.ffn")
        write_fasta({g: s for g, s in proteins.items() if g not in drop}, f"{work}/{label}.faa")
        sample_sheet.append(f"{label},{work}/{label}.ffn,{work}/{label}.faa")
open(f"{work}/samples.csv", "w").write("\n".join(sample_sheet) + "\n")
# the comparison without subsampling needs the complete genomes only
open(f"{work}/samples_all.csv", "w").write("\n".join(sample_sheet[:1] + [line for line in sample_sheet if "_all," in line]) + "\n")


def run_fmh_dnds(sample_sheet_path, scale_factor, out_dir):
    """Run fmh_dnds.py and return its result table (it writes fmh_omega_7.csv into out_dir)."""
    os.makedirs(out_dir, exist_ok=True)
    for old in ("fmh_omega_7.csv", "data.zip"):
        if os.path.exists(f"{out_dir}/{old}"):
            os.remove(f"{out_dir}/{old}")
    command = ["python", "external/dnds-using-fmh/fmh_dnds.py", "--fasta_input_list", sample_sheet_path,
               "--ksize", "7", "--scaled_input", str(scale_factor), "--directory", out_dir,
               "--cores", THREADS, "--threshold", "0"]
    log = subprocess.run(command, capture_output=True, text=True)
    open(f"{out_dir}/fmh_dnds.log", "w").write(log.stdout + log.stderr)
    if log.returncode != 0:
        sys.exit(f"fmh_dnds.py failed; see {out_dir}/fmh_dnds.log")
    return pd.read_csv(f"{out_dir}/fmh_omega_7.csv")


def a_versus_b(results, version):
    """The row of the result table that compares A and B of one version, as a dict."""
    wanted = {f"A_{version}", f"B_{version}"}
    rows = results[[{q, m} == wanted for q, m in zip(results.query_name, results.match_name)]]
    if len(rows) == 0:        # no shared k-mers left: fmh_dnds.py reports nothing for the pair
        return dict(containment_nt=0.0, containment_aa=0.0, fmh_dnds=float("nan"))
    row = rows.iloc[0]
    return dict(containment_nt=row.DNA_max_Cfrac, containment_aa=row.AA_max_Cfrac, fmh_dnds=row["dN/dS"])


# 2. sketches with a scale factor of 500, all versions
sketched = run_fmh_dnds(f"{work}/samples.csv", 500, f"{work}/scaled_500")
rows = [dict(genes=version, **a_versus_b(sketched, version)) for version in VERSIONS]

# 3. every k-mer (scale factor 1), complete genomes only
every_kmer = run_fmh_dnds(f"{work}/samples_all.csv", 1, f"{work}/scaled_1")
rows.append(dict(genes="all_every_kmer", **a_versus_b(every_kmer, "all")))

table = pd.DataFrame(rows)
table.to_csv(f"{pair_dir}/results/fmh_dnds.tsv", sep="\t", index=False)
print(table.to_string(index=False))
