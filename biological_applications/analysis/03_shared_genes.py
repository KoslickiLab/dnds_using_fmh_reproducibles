#!/usr/bin/env python
"""Step 3. Find the genes that are nearly identical between the two genomes of a pair.

    usage:  python analysis/03_shared_genes.py <pair folder>

FracMinHash dN/dS is computed from the k-mers that two genomes share. Here we ask WHERE those
shared k-mers are: for every gene we count how many of its nucleotide 21-mers also occur in the
coding sequences of the other genome. A gene is called "near-identical" when at least half of
its 21-mers are shared (and it has at least 100 of them). Half of the 21-mers surviving
corresponds to roughly 97% nucleotide identity (0.5^(1/21) = 0.968).

Output (in <pair folder>/results/):
    genes_A.tsv, genes_B.tsv   one row per gene: gene, contig, start, end, strand, n_kmers,
                               n_shared, fraction_shared, near_identical
    genome_stats.tsv           one row per genome: size, number of contigs, G+C content, genes
    contigs_A.tsv              contigs of genome A with their lengths (used to draw the figures)
"""
import sys
import os
import pandas as pd
from common import pair_accessions, read_fasta, read_prodigal_genes, gc_percent

K = 21                       # nucleotide k-mer size used by FracMinHash dN/dS (k_aa = 7)
MIN_FRACTION_SHARED = 0.5    # a gene is near-identical if at least this fraction of 21-mers is shared ...
MIN_KMERS = 100              # ... and it has at least this many 21-mers

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def canonical_kmers(sequence):
    """All 21-mers of a sequence. Each k-mer is stored as the alphabetically smaller of itself
    and its reverse complement, so a k-mer and its opposite strand count as the same thing.
    k-mers that contain a letter other than A, C, G and T (for example N) are left out
    (sourmash does both, too)."""
    kmers = []
    sequence = sequence.upper()
    for i in range(len(sequence) - K + 1):
        kmer = sequence[i:i + K]
        if set(kmer) <= set("ACGT"):
            reverse = kmer.translate(COMPLEMENT)[::-1]
            kmers.append(min(kmer, reverse))
    return kmers


pair_dir = sys.argv[1]
os.makedirs(f"{pair_dir}/results", exist_ok=True)
accession_a, accession_b = pair_accessions(pair_dir)
sides = (("A", accession_a), ("B", accession_b))

# k-mers of every gene, for both genomes
gene_kmers = {}
all_kmers = {}
for side, accession in sides:
    coding_sequences = read_fasta(f"data/genes/{accession}.ffn")
    gene_kmers[side] = {gene: canonical_kmers(seq) for gene, seq in coding_sequences.items()}
    all_kmers[side] = set(kmer for kmers in gene_kmers[side].values() for kmer in kmers)

shared = all_kmers["A"] & all_kmers["B"]
print(f"distinct 21-mers: A {len(all_kmers['A']):,}, B {len(all_kmers['B']):,}, shared {len(shared):,}")

stats = []
for side, accession in sides:
    genes = read_prodigal_genes(f"data/genes/{accession}.faa")
    genes["n_kmers"] = [len(gene_kmers[side][g]) for g in genes.gene]
    genes["n_shared"] = [sum(kmer in shared for kmer in gene_kmers[side][g]) for g in genes.gene]
    genes["fraction_shared"] = (genes.n_shared / genes.n_kmers.clip(lower=1)).round(4)
    genes["near_identical"] = (genes.fraction_shared >= MIN_FRACTION_SHARED) & (genes.n_kmers >= MIN_KMERS)
    genes.to_csv(f"{pair_dir}/results/genes_{side}.tsv", sep="\t", index=False)
    near = genes[genes.near_identical]
    print(f"genome {side} ({accession}): {len(genes)} genes, {len(near)} near-identical; "
          f"they hold {near.n_shared.sum() / max(genes.n_shared.sum(), 1):.1%} of the gene k-mers shared with the other genome")

    genome = read_fasta(f"data/genomes/{accession}/genome.fna")
    stats.append(dict(genome=side, accession=accession, size_bp=sum(len(s) for s in genome.values()),
                      n_contigs=len(genome), gc_percent=round(gc_percent("".join(genome.values())), 2), n_genes=len(genes)))
    if side == "A":
        contigs = pd.DataFrame(dict(contig=list(genome), length_bp=[len(s) for s in genome.values()]))
        contigs.sort_values("length_bp", ascending=False).to_csv(f"{pair_dir}/results/contigs_A.tsv", sep="\t", index=False)
pd.DataFrame(stats).to_csv(f"{pair_dir}/results/genome_stats.tsv", sep="\t", index=False)
