#!/usr/bin/env python
"""Step 7. Turn the A-versus-B alignment into a list of shared regions and annotate their genes.

    usage:  python analysis/07_shared_regions.py <pair folder>

1. Keep alignment blocks that are at least 1 kb long and at least 95% identical.
2. Merge blocks that lie within 10 kb of each other on the same contig of genome A into one
   region. Regions are named E1, E2, ... in genome order (contigs from longest to shortest).
3. List the genes NCBI annotated inside each region and sort their products into a few broad
   functional classes by keyword (see FUNCTION_KEYWORDS). The classes only colour the figure
   and support the description in the text; they are not a curated annotation.
   (Genomes without an NCBI annotation get an empty gene list.)
4. For each contig of genome B that carries shared DNA, report how much of it is NOT shared.
   Shared DNA sitting inside a long contig with other sequence around it is integrated in that
   genome; a contig that is shared from end to end could also be a contaminant.

Output (in <pair folder>/results/):
    blocks.tsv          the filtered alignment blocks, with the region each belongs to
    regions.tsv         one row per region: position in A, size, aligned bp, identity, G+C content,
                        length of the contig that carries it, gene counts
    region_genes.tsv    one row per annotated gene in a region, with its functional class
    partner_contigs.tsv one row per contig of genome B that carries shared DNA
"""
import sys
import os
import re
import pandas as pd
from common import pair_accessions, read_coords, read_gff, read_fasta, gc_percent, covered_bases, work_dir

MIN_BLOCK_LENGTH = 1000     # bp
MIN_BLOCK_IDENTITY = 95.0   # percent
MAX_GAP = 10000             # bp; blocks closer than this are merged into one region

# Functional classes, tested in this order; the first class with a matching keyword wins.
FUNCTION_KEYWORDS = {
    "resistance": r"tetracycline|[Tt]et\(|aminoglycoside|lincosamide|macrolide|[Ee]rm\(|rRNA adenine N-6-methyltransferase|"
                  r"beta-lactamase|chloramphenicol|streptothricin|vancomycin|bleomycin|antibiotic|multidrug|"
                  r"mercur|arsen|copper|cadmium|cobalt|heavy metal|chromate|tellur|silver|resistance",
    "conjugation": r"relaxase|[Mm]obiliz|[Mm]obilis|conjug|type IV secret|type-IV secret|VirB|VirD|Tra[A-Z]\b|TrbL|"
                   r"T4SS|coupling protein|Mob[ACQV]\b|PcfJ|PcfB|TcpE|nickase|plasmid transfer|PrgI|Maff2|"
                   r"mobile element protein|transposon-transfer|TnpV|TnpW|TnpX",
    "integration": r"integrase|recombinase|excisionase|resolvase|site-specific",
    "transposase": r"transposase|IS\d+ family|insertion sequence",
    "defence": r"restriction|BREX|Brx|Pgl[XZ]|toxin|antitoxin|abortive|CRISPR|DNA[- ]methyltransferase|DNA methylase|"
               r"DNA \(cytosine|adenine-specific|phosphorothioat|anti-?restriction|ArdA|DarA|addiction",
    "phage": r"phage|capsid|terminase|portal|tail|holin|prohead",
    "replication": r"replication initiat|Rep[ABC]\b|Par[AB]\b|partition|plasmid|primase|helicase",
    "hypothetical": r"hypothetical|DUF\d+|uncharacteri[sz]ed|unknown function",
}


def functional_class(product):
    for name, keywords in FUNCTION_KEYWORDS.items():
        if re.search(keywords, product):
            return name
    return "other"


pair_dir = sys.argv[1]
accession_a, accession_b = pair_accessions(pair_dir)
results = f"{pair_dir}/results"
genome_a = read_fasta(f"data/genomes/{accession_a}/genome.fna")

# 1. filter the alignment blocks (in genome A, start_ref < end_ref always)
blocks = read_coords(f"{work_dir(pair_dir, 'alignments')}/{accession_b}.coords")
blocks = blocks[(blocks.len_ref >= MIN_BLOCK_LENGTH) & (blocks.len_qry >= MIN_BLOCK_LENGTH)
                & (blocks.identity >= MIN_BLOCK_IDENTITY)]
# order: contigs of genome A from longest to shortest, then by position (the order in which the figure shows them)
blocks = blocks.sort_values(["contig_len_ref", "contig_ref", "start_ref"], ascending=[False, True, True]).reset_index(drop=True)

# 2. merge neighbouring blocks into regions
regions = []               # each region: dict(contig, start, end)
region_of_block = []
for block in blocks.itertuples():
    last = regions[-1] if regions else None
    if last and last["contig"] == block.contig_ref and block.start_ref <= last["end"] + MAX_GAP:
        last["end"] = max(last["end"], block.end_ref)
    else:
        regions.append(dict(contig=block.contig_ref, start=block.start_ref, end=block.end_ref))
    region_of_block.append(f"E{len(regions)}")
blocks["region"] = region_of_block
blocks.to_csv(f"{results}/blocks.tsv", sep="\t", index=False)

# 3. genes of genome A inside each region
gff = f"data/genomes/{accession_a}/annotation.gff"
if os.path.exists(gff):
    annotation = read_gff(gff)
else:
    annotation = pd.DataFrame(columns=["contig", "start", "end", "strand", "type", "locus_tag", "product"])
annotation["function"] = annotation["product"].map(functional_class)
region_rows = []
gene_tables = []
for number, region in enumerate(regions, start=1):
    name = f"E{number}"
    mine = blocks[blocks.region == name]
    aligned = covered_bases(list(zip(mine.start_ref, mine.end_ref)), region["start"], region["end"])
    genes = annotation[(annotation.contig == region["contig"]) & (annotation.end >= region["start"])
                       & (annotation.start <= region["end"])].copy()
    genes.insert(0, "region", name)
    gene_tables.append(genes)
    counts = genes.function.value_counts()
    row = dict(region=name, contig=region["contig"], start=region["start"], end=region["end"],
               span_bp=region["end"] - region["start"] + 1, aligned_bp=aligned,
               contig_length_bp=len(genome_a[region["contig"]]),
               gc_percent=round(gc_percent(genome_a[region["contig"]][region["start"] - 1:region["end"]]), 1),
               identity=round((mine.identity * mine.len_ref).sum() / mine.len_ref.sum(), 3),   # length-weighted mean
               longest_block_bp=int(mine.len_ref.max()),
               partner_contigs=",".join(sorted(set(mine.contig_qry))),
               n_genes=len(genes))
    for function in list(FUNCTION_KEYWORDS) + ["other"]:
        row[f"n_{function}"] = int(counts.get(function, 0))
    region_rows.append(row)
region_table = pd.DataFrame(region_rows)
region_table.to_csv(f"{results}/regions.tsv", sep="\t", index=False)
pd.concat(gene_tables).to_csv(f"{results}/region_genes.tsv", sep="\t", index=False)

# 4. the contigs of genome B that carry the shared DNA
partner_rows = []
for contig, mine in blocks.groupby("contig_qry"):
    length = int(mine.contig_len_qry.iloc[0])
    intervals = [(min(s, e), max(s, e)) for s, e in zip(mine.start_qry, mine.end_qry)]
    shared = covered_bases(intervals, 1, length)
    partner_rows.append(dict(contig=contig, length_bp=length, shared_bp=shared, not_shared_bp=length - shared,
                             regions=",".join(sorted(set(mine.region), key=lambda r: int(r[1:])))))
partner_table = pd.DataFrame(partner_rows).sort_values("length_bp", ascending=False)
partner_table.to_csv(f"{results}/partner_contigs.tsv", sep="\t", index=False)

pd.set_option("display.width", 220)
print(region_table[["region", "contig", "start", "end", "span_bp", "aligned_bp", "identity", "n_genes", "n_conjugation",
                    "n_integration", "n_resistance", "n_defence", "n_transposase"]].to_string(index=False))
print(partner_table.head(12).to_string(index=False))
