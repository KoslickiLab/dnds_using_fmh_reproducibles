#!/usr/bin/env python
"""Step 2. Predict protein-coding genes with Prodigal, as GTDB does.

    usage:  python analysis/02_call_genes.py <pair folder>

GTDB (release 214) calls genes with Prodigal 2.6.3 (options -m -p single). The translation
table is the one that NCBI states for the genome. For genomes without an NCBI annotation, genes
are called with tables 11 and 4, and the table-4 genes are kept if their coding density exceeds
that of table 11 by more than 0.05 and is above 0.7. The same is done here, so that the gene
calls equal those from which the FracMinHash dN/dS values of the manuscript were computed
(tested by checks/check_against_gtdb.py). Of the genomes used here, one is called with
table 4 (GCA_923184405.1).

For each accession in the pair's genomes.tsv this creates
    data/genes/<accession>.ffn   coding sequences (nucleotide)
    data/genes/<accession>.faa   proteins
and adds a line to data/genes/translation_tables.tsv (how the table was chosen).
Genomes that already have gene calls are skipped.
"""
import sys
import os
import re
import subprocess
from multiprocessing import Pool
from common import read_genomes, read_fasta, read_prodigal_genes, covered_bases, THREADS


def ncbi_translation_table(accession):
    """The translation table of NCBI's annotation of a genome (transl_table=... in its GFF file), or None."""
    gff = f"data/genomes/{accession}/annotation.gff"
    if os.path.exists(gff):
        for line in open(gff):
            found = re.search(r"transl_table=(\d+)", line)
            if found:
                return int(found.group(1))
    return None


def run_prodigal(accession, table):
    """Call genes with one translation table; returns the coding density (fraction of the genome
    covered by genes, overlapping genes counted once)."""
    genome = f"data/genomes/{accession}/genome.fna"
    prefix = f"data/genes/{accession}.table{table}"
    subprocess.run(["prodigal", "-m", "-p", "single", "-g", str(table), "-q", "-i", genome,
                    "-d", f"{prefix}.ffn", "-a", f"{prefix}.faa", "-o", os.devnull], check=True)
    genes = read_prodigal_genes(f"{prefix}.faa")
    coding = sum(covered_bases(list(zip(on_contig.start, on_contig.end)), 1, on_contig.end.max())
                 for _, on_contig in genes.groupby("contig"))
    return coding / sum(len(sequence) for sequence in read_fasta(genome).values())


def call_genes(accession):
    stated = ncbi_translation_table(accession)
    if stated is not None:
        tables, chosen, reason = [stated], stated, "stated by NCBI"
        density = {stated: run_prodigal(accession, stated)}
    else:
        tables = [11, 4]
        density = {table: run_prodigal(accession, table) for table in tables}
        chosen = 4 if density[4] - density[11] > 0.05 and density[4] > 0.7 else 11
        reason = f"coding density {density[11]:.3f} with table 11, {density[4]:.3f} with table 4"
    for table in tables:
        for extension in ("ffn", "faa"):
            if table == chosen:
                os.replace(f"data/genes/{accession}.table{table}.{extension}", f"data/genes/{accession}.{extension}")
            else:
                os.remove(f"data/genes/{accession}.table{table}.{extension}")
    return f"{accession}\t{chosen}\t{reason}\n"


if __name__ == "__main__":
    pair_dir = sys.argv[1]
    os.makedirs("data/genes", exist_ok=True)
    accessions = sorted(set(read_genomes(pair_dir).accession))
    todo = [accession for accession in accessions if not os.path.exists(f"data/genes/{accession}.faa")]
    print(f"calling genes in {len(todo)} genomes ({len(accessions) - len(todo)} already done)")
    with Pool(int(THREADS)) as pool:
        lines = pool.map(call_genes, todo)
    new_file = not os.path.exists("data/genes/translation_tables.tsv")
    with open("data/genes/translation_tables.tsv", "a") as out:
        if new_file:
            out.write("accession\ttranslation_table\tchosen_because\n")
        out.writelines(lines)
