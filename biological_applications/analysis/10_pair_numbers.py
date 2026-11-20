#!/usr/bin/env python
"""Step 10. Collect every number that the manuscript quotes about a pair into LaTeX macros.

    usage:  python analysis/10_pair_numbers.py <pair folder>

Reads the result tables of steps 3-9 (no genome data needed) and writes
output/numbers/<pair folder name>.tex, one \\providecommand per number. The manuscript \\input's
that file, so the text cannot disagree with the analysis. Macro names start with the pair's
prefix (settings.ini), e.g. \\FRFmhDnds.
"""
import sys
import os
import pandas as pd
from common import read_genomes, read_settings, write_numbers, thousands, kb

NUMBER_WORDS = ["One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Ten", "Eleven", "Twelve",
                "Thirteen", "Fourteen", "Fifteen", "Sixteen", "Seventeen", "Eighteen", "Nineteen", "Twenty",
                "TwentyOne", "TwentyTwo", "TwentyThree", "TwentyFour", "TwentyFive", "TwentySix", "TwentySeven",
                "TwentyEight", "TwentyNine", "Thirty"]
PRESENT = 0.9      # a region counts as present in a control genome if at least 90% of it is covered

pair_dir = os.path.normpath(sys.argv[1])
results = f"{pair_dir}/results"
settings = read_settings(pair_dir)
genomes = read_genomes(pair_dir)

fmh = pd.read_csv(f"{results}/fmh_dnds.tsv", sep="\t").set_index("genes")
genes = {side: pd.read_csv(f"{results}/genes_{side}.tsv", sep="\t") for side in "AB"}
orthologs = pd.read_csv(f"{results}/orthologs.tsv", sep="\t")
regions = pd.read_csv(f"{results}/regions.tsv", sep="\t")
blocks = pd.read_csv(f"{results}/blocks.tsv", sep="\t")
partners = pd.read_csv(f"{results}/partner_contigs.tsv", sep="\t")
presence = pd.read_csv(f"{results}/presence.tsv", sep="\t")
genome_stats = pd.read_csv(f"{results}/genome_stats.tsv", sep="\t").set_index("genome")


numbers = []      # (macro name without prefix, value, meaning)


def two_decimals(value):
    """A dN/dS value as text; pairs that share no k-mers after gene removal have no estimate."""
    return "no estimate" if pd.isna(value) else f"{value:.2f}"


# ---- the FracMinHash estimate, with and without the near-identical genes
numbers += [
    ("FmhDnds", two_decimals(fmh.loc['all', 'fmh_dnds']), "FracMinHash dN/dS of the pair, all genes"),
    ("FmhDndsWithout", two_decimals(fmh.loc['without_near_identical', 'fmh_dnds']), "FracMinHash dN/dS after removing near-identical genes from both genomes"),
    ("FmhDndsWithoutLoose", two_decimals(fmh.loc['without_genes_sharing_25_percent', 'fmh_dnds']), "the same when genes sharing at least 25 percent of their 21-mers are removed"),
    ("FmhDndsWithoutStrict", two_decimals(fmh.loc['without_genes_sharing_75_percent', 'fmh_dnds']), "the same when only genes sharing at least 75 percent of their 21-mers are removed"),
    ("FmhDndsEveryKmer", two_decimals(fmh.loc['all_every_kmer', 'fmh_dnds']), "FracMinHash dN/dS of the pair computed from every k-mer (scale factor 1)"),
    ("ContainmentNt", f"{fmh.loc['all', 'containment_nt']:.3f}", "nucleotide 21-mer containment, all genes"),
    ("ContainmentAa", f"{fmh.loc['all', 'containment_aa']:.3f}", "protein 7-mer containment, all genes"),
    ("ContainmentNtWithout", f"{fmh.loc['without_near_identical', 'containment_nt']:.3f}", "nucleotide containment without near-identical genes"),
    ("ContainmentAaWithout", f"{fmh.loc['without_near_identical', 'containment_aa']:.3f}", "protein containment without near-identical genes"),
]

# ---- the two genomes and their near-identical genes
for side in "AB":
    table = genes[side]
    near = table[table.near_identical]
    stats = genome_stats.loc[side]
    numbers += [
        (f"GenomeMb{side}", f"{stats.size_bp / 1e6:.2f}", f"genome {side}: assembly size (Mb)"),
        (f"Contigs{side}", str(stats.n_contigs), f"genome {side}: number of contigs or scaffolds"),
        (f"Gc{side}", f"{stats.gc_percent:.1f}", f"genome {side}: G+C content (percent)"),
        (f"Genes{side}", thousands(len(table)), f"genome {side}: predicted genes"),
        (f"NearGenes{side}", str(len(near)), f"genome {side}: near-identical genes (>= 50 percent of 21-mers shared)"),
        (f"NearPercentGenes{side}", f"{100 * len(near) / len(table):.1f}", f"genome {side}: near-identical genes as percent of all genes"),
        (f"NearKb{side}", kb((near.end - near.start + 1).sum()), f"genome {side}: total length of near-identical genes (kb)"),
        (f"NearShare{side}", f"{100 * near.n_shared.sum() / table.n_shared.sum():.0f}",
         f"genome {side}: percent of all shared gene 21-mers that lie in near-identical genes"),
    ]

# ---- CheckM estimates of the two genomes, for genomes that are GTDB R214 species representatives
gtdb = pd.read_csv("all_pairs/genomes.tsv", sep="\t").set_index("accession")
for side in "AB":
    accession = genome_stats.loc[side, "accession"]
    if accession in gtdb.index:
        numbers += [(f"CheckmCompleteness{side}", f"{gtdb.checkm_completeness[accession]:.1f}", f"genome {side}: CheckM completeness in the GTDB R214 metadata (percent)"),
                    (f"CheckmContamination{side}", f"{gtdb.checkm_contamination[accession]:.2f}", f"genome {side}: CheckM contamination in the GTDB R214 metadata (percent)")]

# ---- conventional dN/dS of the orthologs
background = orthologs[~orthologs.near_identical]
estimable = background[background.dNdS.notna() & (background.dS > 0)]
numbers += [
    ("Orthologs", thousands(len(orthologs)), "reciprocal-best-hit orthologs that could be aligned"),
    ("NearOrthologs", str(int(orthologs.near_identical.sum())), "orthologs that are near-identical genes"),
    ("BackgroundOrthologs", thousands(len(background)), "orthologs that are not near-identical"),
    ("BackgroundEstimable", thousands(len(estimable)), "background orthologs with an estimable dN/dS"),
    ("AaIdentity", f"{background.aa_identity.median():.0f}", "median amino-acid identity of background orthologs (percent)"),
    ("ConvDnds", f"{estimable.dNdS.median():.2f}", "conventional dN/dS (YN), median over background orthologs"),
    ("ConvDs", f"{estimable.dS.median():.1f}", "conventional dS (YN), median over background orthologs"),
    ("ConvAboveOne", str(int((estimable.dNdS > 1).sum())), "background orthologs with conventional dN/dS above 1"),
]

# ---- the shared regions
genome_a_bp = genome_stats.loc["A", "size_bp"]
longest = blocks.loc[blocks.len_ref.idxmax()]
numbers += [
    ("Regions", str(len(regions)), "number of shared regions (blocks >= 1 kb and >= 95 percent identical, merged within 10 kb)"),
    ("RegionsKb", kb(regions.span_bp.sum()), "total span of the shared regions in genome A (kb)"),
    ("RegionsPercentA", f"{100 * regions.span_bp.sum() / genome_a_bp:.1f}", "shared regions as percent of genome A"),
    ("AlignedKb", kb(regions.aligned_bp.sum()), "aligned at >= 95 percent identity (kb)"),
    ("AlignedNearPerfectKb", kb(blocks[blocks.identity >= 99.9].len_ref.sum()), "aligned at >= 99.9 percent identity (kb)"),
    ("LongestBlockKb", kb(longest.len_ref, 1), "longest single alignment block (kb)"),
    ("LongestBlockIdentity", f"{min(longest.identity, 99.99):.2f}", "identity of the longest block (percent; capped at 99.99)"),
    ("PartnerContigs", str(len(partners)), "contigs of genome B that carry shared DNA"),
    ("PartnerContigsEmbedded", str(int((partners.not_shared_bp >= 10000).sum())), "of those, contigs with >= 10 kb of other sequence as well"),
    ("LargestPartnerContigKb", kb(partners.length_bp.max()), "largest contig of genome B that carries shared DNA (kb)"),
    ("PartnerContigsWhollyShared", str(int((partners.not_shared_bp == 0).sum())), "contigs of genome B that consist of shared DNA only"),
]
# scaffold gaps (runs of N) inside the regions or within 20 kb of them, in each genome (step 8)
gaps = pd.read_csv(f"{results}/gaps_near_regions.tsv", sep="\t")
for side in "AB":
    mine = gaps[gaps.genome == side]
    n_gaps = int(mine[["gaps_left_20kb", "gaps_inside", "gaps_right_20kb"]].sum().sum())
    numbers.append((f"GapsNear{side}", str(n_gaps), f"genome {side}: scaffold gaps inside or within 20 kb of the shared regions"))
numbers += [("RegionsGcMin", f"{regions.gc_percent.min():.1f}", "lowest G+C content among the shared regions (percent)"),
            ("RegionsGcMax", f"{regions.gc_percent.max():.1f}", "highest G+C content among the shared regions (percent)"),
            ("RegionsIdentityMin", f"{regions.identity.min():.1f}", "lowest identity to genome B among the shared regions (percent)"),
            ("RegionsIdentityMax", f"{min(regions.identity.max(), 99.9):.1f}", "highest identity to genome B among the shared regions (percent; capped at 99.9)")]
for word, region in zip(NUMBER_WORDS, regions.itertuples()):
    numbers += [
        (f"Region{word}Kb", kb(region.span_bp), f"{region.region}: span (kb)"),
        (f"Region{word}Identity", f"{min(region.identity, 99.99):.2f}", f"{region.region}: identity to genome B (percent; capped at 99.99)"),
        (f"Region{word}Gc", f"{region.gc_percent:.1f}", f"{region.region}: G+C content (percent)"),
        (f"Region{word}Genes", str(region.n_genes), f"{region.region}: annotated genes"),
        (f"Region{word}Integrases", str(region.n_integration), f"{region.region}: genes annotated as integrase or recombinase"),
        (f"Region{word}ContigKb", thousands(region.contig_length_bp / 1000), f"{region.region}: length of the genome-A contig or scaffold that carries it (kb)"),
    ]

# ---- presence of the regions in genome B and the control genomes
for genome in genomes[genomes.role != "A"].itertuples():
    mine = presence[presence.accession == genome.accession]
    present = mine[mine.covered >= PRESENT]
    numbers.append((f"Present{genome.key}", str(len(present)),
                    f"regions at least 90 percent present in {genome.label} ({genome.accession})"))
    if len(present):
        numbers.append((f"SmallestHostContigKb{genome.key}", thousands(present.located_on_length.min() / 1000),
                        f"smallest contig of {genome.label} that carries one of these regions (kb)"))
# single cells of the presence table that the text quotes; settings.ini:  quoted_cells = Key:E3, Key:E5
for cell in settings.get("quoted_cells", "").split(","):
    if ":" not in cell:
        continue
    key, region_name = cell.strip().split(":")
    accession = genomes[genomes.key == key].accession.iloc[0]
    row = presence[(presence.accession == accession) & (presence.region == region_name)].iloc[0]
    word = NUMBER_WORDS[int(region_name[1:]) - 1]
    numbers += [(f"Covered{key}Region{word}", f"{100 * row.covered:.0f}", f"{region_name}: percent present in {key}"),
                (f"Identity{key}Region{word}", f"{min(row.identity, 99.99):.2f}", f"{region_name}: identity in {key} (percent)"),
                (f"Flank{key}Region{word}", f"{100 * row.flank_covered:.0f}", f"{region_name}: percent of the 10-kb flanks aligned in {key}")]

# counts that the text quotes; settings.ini:  quoted_counts = Name=E1:KeyOne+KeyTwo; OtherName=E2:KeyThree
# -> \<prefix>Count<Name> = how many of the listed genomes carry at least 90 percent of the region,
#    \<prefix>Count<Name>Of = how many genomes were listed
for item in settings.get("quoted_counts", "").split(";"):
    if "=" not in item:
        continue
    name, rest = item.strip().split("=")
    region_name, keys = rest.split(":")
    accessions = [genomes[genomes.key == key.strip()].accession.iloc[0] for key in keys.split("+")]
    mine = presence[(presence.region == region_name.strip()) & presence.accession.isin(accessions)]
    numbers += [(f"Count{name}", str(int((mine.covered >= PRESENT).sum())), f"{region_name}: genomes among {keys} with at least 90 percent of the region"),
                (f"Count{name}Of", str(len(accessions)), f"{region_name}: number of genomes in that list")]

# ---- read depth and read pairs at the shared region (optional steps 16 and 17; only for the pair that has them)
if os.path.exists(f"{results}/read_depth.tsv"):
    depth = pd.read_csv(f"{results}/read_depth.tsv", sep="\t").set_index("sequence").median_depth
    pairs = pd.read_csv(f"{results}/read_pairs.tsv", sep="\t").set_index(["inside_read", "other_read_lies_on"]).read_pairs
    at_end, interior = pairs["at an end of the region"], pairs["interior of the region"]
    numbers += [
        ("DepthRegion", f"{depth['shared region']:.0f}", "median read depth of the shared region in the source metagenome of genome B"),
        ("DepthFlank", f"{depth['rest of the scaffold that carries the region (genome A)']:.0f}", "median read depth of the rest of the genome-A scaffold that carries the region"),
        ("DepthGenomeB", f"{depth['genome B (contigs that consist of shared sequence left out)']:.0f}", "median read depth of genome B"),
        ("DepthExpectedBoth", f"{depth['expected for the region if both genomes carried it']:.0f}", "depth expected for the region if both genomes carried it"),
        ("PairsEnd", str(int(at_end.sum())), "read pairs with one read within 500 bp of an end of the region and the other read outside the region"),
        ("PairsEndFlank", str(int(at_end["flanking sequence of the same scaffold (genome A)"])), "of those, pairs whose other read lies in the flanking sequence of genome A"),
        ("PairsEndGenomeB", str(int(at_end["a contig of genome B"])), "of those, pairs whose other read lies on a contig of genome B"),
        ("PairsInteriorGenomeB", str(int(interior["a contig of genome B"])), "read pairs with one read in the interior of the region and the other on a contig of genome B (chimeric pairs)"),
    ]

write_numbers(f"output/numbers/{os.path.basename(pair_dir)}.tex", settings["prefix"], numbers)
