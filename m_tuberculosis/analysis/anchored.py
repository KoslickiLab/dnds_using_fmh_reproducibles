"""Reference-anchored coding sequences of the isolates: the input construction, and nothing else.

Every isolate is represented the way Gan et al. 2025 represent it: the inferred MTBC ancestral
sequence (same coordinates as H37Rv, NC_000962.3) with that isolate's substitutions applied. From it
the coding sequences of the reference annotation are cut out and translated (table 11). No k-mers,
hashes or containments are computed here; that is done by the FracMinHash dN/dS tool
(external/dnds-using-fmh) from the FASTA files these functions produce.

    rows of results/groups/gene_groups.tsv      one per coding sequence (or high-quality segment)
    work/variants/<sample>.tsv.gz               substitutions found from an assembly (step 14)
    work/snv_reads/<run>.calls.tsv.gz           substitutions found from reads (steps 07-08)
"""
import gzip
import os
import pandas as pd

MTB = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GROUPS = ["all_cds", "essential", "nonessential", "pe_ppe", "pe_ppe_hq"]
COMPLEMENT = str.maketrans("ACGT", "TGCA")
CODONS = {}
for a, first in enumerate("TCAG"):
    for b, second in enumerate("TCAG"):
        for c, third in enumerate("TCAG"):
            CODONS[first + second + third] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[a * 16 + b * 4 + c]


def translate(dna):
    """Translation table 11 (identical to the standard table over the 64 codons); a trailing stop is dropped."""
    protein = "".join([CODONS.get(dna[i:i + 3], "X") for i in range(0, len(dna) - len(dna) % 3, 3)])
    return protein[:-1] if protein.endswith("*") else protein


def read_backbone():
    parts = [line.strip() for line in open(f"{MTB}/data/reference/tb.ancestor.fasta") if not line.startswith(">")]
    return "".join(parts).upper()


def gene_table():
    """Coding sequences (and high-quality PE/PPE segments) whose length is a whole number of codons."""
    genes = pd.read_csv(f"{MTB}/results/groups/gene_groups.tsv", sep="\t")
    return genes[(genes.end - genes.start + 1) % 3 == 0].sort_values("start").reset_index(drop=True)


def substitutions_of(sample, source):
    """{genome position (1-based): alternative base} of one isolate ("assembly" or "reads")."""
    path = (f"{MTB}/work/variants/{sample}.tsv.gz" if source == "assembly"
            else f"{MTB}/work/snv_reads/{sample}.calls.tsv.gz")
    out = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            field = line.split("\t")
            alternative = field[2].strip().upper()
            if len(alternative) == 1 and alternative in "ACGT":
                out[int(field[0])] = alternative
    return out


def apply(backbone, substitutions):
    sequence = list(backbone)
    for position, alternative in substitutions.items():
        if 1 <= position <= len(sequence):
            sequence[position - 1] = alternative
    return "".join(sequence)


def coding_sequences(genome, genes):
    """[(locus_tag, group, in_all_cds, cds, protein)] for one genome, in annotation order."""
    out = []
    for locus_tag, start, end, strand, group, in_all in zip(genes.locus_tag, genes.start, genes.end,
                                                            genes.strand, genes.group, genes.in_all_cds):
        cds = genome[start - 1:end]
        if strand == "-":
            cds = cds.translate(COMPLEMENT)[::-1]
        out.append((locus_tag, group, bool(in_all), cds, translate(cds)))
    return out


def members(group, in_all, group_name):
    """Whether a coding sequence belongs to the named group of sequences."""
    return in_all if group_name == "all_cds" else group == group_name
