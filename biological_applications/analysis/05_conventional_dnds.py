#!/usr/bin/env python
"""Step 5. Conventional, alignment-based dN/dS of the orthologous genes of the pair.

    usage:  python analysis/05_conventional_dnds.py <pair folder>

This is the yardstick the FracMinHash value is compared with.
1. Orthologs = reciprocal best DIAMOND blastp hits between the two proteomes
   (e-value <= 1e-5, alignment covering >= 60% of both proteins).
2. Each ortholog pair is aligned as proteins (global alignment, BLOSUM62, gap open 10,
   gap extend 0.5), and the alignment is copied onto the coding sequences codon by codon.
   Columns with a gap or a stop codon are dropped.
3. KaKs_Calculator (Yang-Nielsen model, genetic code 11) estimates dN and dS for each pair.

Output: <pair folder>/results/orthologs.tsv
    gene_A, gene_B, aa_identity (percent, over aligned columns), n_codons, dN, dS, dNdS,
    near_identical (True if step 3 flagged the genome-A gene)
"""
import sys
import subprocess
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
from common import pair_accessions, read_fasta, work_dir, THREADS

MIN_COVERAGE = 0.6
MAX_EVALUE = 1e-5
MIN_CODONS = 50
STOP_CODONS = {"TAA", "TAG", "TGA"}

pair_dir = sys.argv[1]
work = work_dir(pair_dir, "dnds")
accession_a, accession_b = pair_accessions(pair_dir)
proteins = {"A": f"data/genes/{accession_a}.faa", "B": f"data/genes/{accession_b}.faa"}

# ---- 1. reciprocal best hits
columns = ["query", "subject", "identity", "length", "query_length", "subject_length", "evalue", "bitscore"]


def best_hits(query_side, subject_side):
    """Best hit of every query protein in the other proteome (after the e-value and coverage filters)."""
    out = f"{work}/{query_side}_vs_{subject_side}.tsv"
    subprocess.run(["diamond", "makedb", "--in", proteins[subject_side], "-d", f"{work}/{subject_side}", "--quiet"], check=True)
    subprocess.run(["diamond", "blastp", "--more-sensitive", "--quiet", "--threads", THREADS,
                    "-q", proteins[query_side], "-d", f"{work}/{subject_side}", "-o", out,
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"], check=True)
    hits = pd.read_csv(out, sep="\t", header=None, names=columns)
    hits = hits[(hits.evalue <= MAX_EVALUE) & (hits.length >= MIN_COVERAGE * hits.query_length)
                & (hits.length >= MIN_COVERAGE * hits.subject_length)]
    hits = hits.sort_values("bitscore", ascending=False).drop_duplicates("query")
    return dict(zip(hits["query"], hits["subject"]))


a_to_b = best_hits("A", "B")
b_to_a = best_hits("B", "A")
orthologs = [(a, b) for a, b in a_to_b.items() if b_to_a.get(b) == a]
print(f"{len(orthologs)} reciprocal best hits")

# ---- 2. protein alignment -> codon alignment
aligner = PairwiseAligner()
aligner.mode = "global"
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

protein_a, protein_b = read_fasta(proteins["A"]), read_fasta(proteins["B"])
cds_a, cds_b = read_fasta(f"data/genes/{accession_a}.ffn"), read_fasta(f"data/genes/{accession_b}.ffn")

rows = []
with open(f"{work}/orthologs.axt", "w") as axt:
    for gene_a, gene_b in orthologs:
        seq_a, seq_b = protein_a[gene_a].rstrip("*"), protein_b[gene_b].rstrip("*")
        alignment = aligner.align(seq_a, seq_b)[0]
        codons_a, codons_b, identical = [], [], 0
        # alignment.aligned lists the gap-free stretches as ((start_a, end_a), (start_b, end_b))
        for (start_a, end_a), (start_b, end_b) in zip(*alignment.aligned):
            for offset in range(end_a - start_a):
                i, j = start_a + offset, start_b + offset
                codon_a = cds_a[gene_a][3 * i:3 * i + 3].upper()
                codon_b = cds_b[gene_b][3 * j:3 * j + 3].upper()
                if codon_a in STOP_CODONS or codon_b in STOP_CODONS or set(codon_a + codon_b) - set("ACGT"):
                    continue
                codons_a.append(codon_a)
                codons_b.append(codon_b)
                identical += seq_a[i] == seq_b[j]
        if len(codons_a) < MIN_CODONS:
            continue
        name = f"{gene_a}__{gene_b}"
        axt.write(f"{name}\n{''.join(codons_a)}\n{''.join(codons_b)}\n\n")
        rows.append(dict(name=name, gene_A=gene_a, gene_B=gene_b, n_codons=len(codons_a),
                         aa_identity=round(100 * identical / len(codons_a), 2)))
table = pd.DataFrame(rows)

# ---- 3. dN and dS with KaKs_Calculator
subprocess.run(["KaKs_Calculator", "-i", f"{work}/orthologs.axt", "-o", f"{work}/orthologs.kaks", "-c", "11", "-m", "YN"],
               stdout=subprocess.DEVNULL, check=True)
kaks = pd.read_csv(f"{work}/orthologs.kaks", sep="\t")[["Sequence", "Ka", "Ks", "Ka/Ks"]]
kaks.columns = ["name", "dN", "dS", "dNdS"]
for column in ("dN", "dS", "dNdS"):
    kaks[column] = pd.to_numeric(kaks[column], errors="coerce")     # KaKs_Calculator writes "NA" when not estimable
table = table.merge(kaks, on="name", how="left").drop(columns="name")

near_identical = set(pd.read_csv(f"{pair_dir}/results/genes_A.tsv", sep="\t").query("near_identical").gene)
table["near_identical"] = table.gene_A.isin(near_identical)
table.to_csv(f"{pair_dir}/results/orthologs.tsv", sep="\t", index=False)

estimable = table[table.dNdS.notna() & (table.dS > 0)]
background = estimable[~estimable.near_identical]
print(f"{len(table)} aligned orthologs, {len(estimable)} with an estimable dN/dS")
print(f"background orthologs: median amino-acid identity {table[~table.near_identical].aa_identity.median():.1f}%, "
      f"median dN/dS {background.dNdS.median():.3f}, median dS {background.dS.median():.2f}")
print(f"near-identical orthologs: {int(table.near_identical.sum())}")
