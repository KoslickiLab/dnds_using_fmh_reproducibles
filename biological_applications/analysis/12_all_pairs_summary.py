#!/usr/bin/env python
"""Step 12. One row per genome pair of the hierarchical edge bundling figure: where its shared
k-mers lie, what remains of its FracMinHash dN/dS without them, and what kind of pair it is.

    usage:  python analysis/12_all_pairs_summary.py          (after analysis/11_all_pairs_run.sh)

Inputs
    all_pairs/pairs.tsv      the 140 pairs and their FracMinHash dN/dS in the manuscript's figure
    all_pairs/genomes.tsv    GTDB R214 metadata of the 243 genomes (taxonomy, isolate / MAG / SAG,
                             CheckM completeness and contamination)
    work/all_pairs/<pair>/results/   genes_A.tsv, genes_B.tsv, fmh_dnds.tsv, orthologs.tsv (steps 3-5)

Classes of pairs (column "class")
    A contig is called "wholly shared" when at least 80% of its genes are near-identical to the
    other genome of the pair: the contig as a whole is DNA of the other genome.
      related genomes
            fewer than half of the pair's shared 21-mers lie in near-identical genes: the
            similarity is spread over the genomes
      otherwise, if a metagenome-assembled genome (MAG) or single-cell genome (SAG) is involved:
        MAG or SAG, foreign contigs
            at least 20% of the near-identical genes of a MAG/SAG lie on wholly shared contigs,
            and fewer than half of the orthologs are near-identical (diverged background)
        MAG or SAG, close relatives with partial overlap
            the same, but at least half of the orthologs are near-identical
        MAG or SAG, unresolved
            the shared DNA lies inside longer contigs
      otherwise (both genomes are isolates):
        isolates with wholly shared contigs
            at least half of the near-identical genes of either genome lie on wholly shared contigs
        isolates with embedded shared DNA
            the shared DNA lies inside longer contigs

Output: all_pairs/results/all_pairs_summary.tsv
"""
import pandas as pd

WHOLLY_SHARED = 0.8       # a contig is wholly shared if at least this fraction of its genes is near-identical
RANKS = ["domain", "phylum", "class", "order", "family", "genus"]

pairs = pd.read_csv("all_pairs/pairs.tsv", sep="\t")
genomes = pd.read_csv("all_pairs/genomes.tsv", sep="\t").set_index("accession")


def lowest_shared_rank(taxonomy_1, taxonomy_2):
    """Lowest rank at which two GTDB taxonomy strings (d__...;p__...;...;s__...) agree; 'none' if the domains differ."""
    shared = "none"
    for rank, taxon_1, taxon_2 in zip(RANKS, taxonomy_1.split(";"), taxonomy_2.split(";")):
        if taxon_1 != taxon_2:
            break
        shared = rank
    return shared


rows = []
for pair in pairs.itertuples():
    results = f"work/all_pairs/{pair.pair}/results"
    fmh = pd.read_csv(f"{results}/fmh_dnds.tsv", sep="\t").set_index("genes")
    orthologs = pd.read_csv(f"{results}/orthologs.tsv", sep="\t")
    estimable = orthologs[orthologs.dNdS.notna() & (orthologs.dS > 0)]
    row = dict(pair=pair.pair, genome_1=pair.genome_1, genome_2=pair.genome_2,
               species_1=genomes.gtdb_species[pair.genome_1], species_2=genomes.gtdb_species[pair.genome_2],
               lowest_shared_rank=lowest_shared_rank(genomes.gtdb_taxonomy[pair.genome_1], genomes.gtdb_taxonomy[pair.genome_2]),
               type_1=genomes.genome_type[pair.genome_1], type_2=genomes.genome_type[pair.genome_2],
               containment_nt=fmh.loc["all", "containment_nt"], containment_aa=fmh.loc["all", "containment_aa"],
               fmh_dnds=fmh.loc["all", "fmh_dnds"], fmh_dnds_in_paper=pair.fmh_dnds_in_paper,
               fmh_dnds_every_kmer=fmh.loc["all_every_kmer", "fmh_dnds"],
               fmh_dnds_without_near_identical=fmh.loc["without_near_identical", "fmh_dnds"],
               n_orthologs=len(orthologs), n_near_identical_orthologs=int(orthologs.near_identical.sum()),
               background_aa_identity=orthologs[~orthologs.near_identical].aa_identity.median(),
               conventional_median_dnds=estimable.dNdS.median())

    # where the near-identical genes are, in each of the two genomes
    for side, number in (("A", 1), ("B", 2)):
        genes = pd.read_csv(f"{results}/genes_{side}.tsv", sep="\t")
        near = genes[genes.near_identical]
        # share of this genome's shared 21-mers that lies in near-identical genes
        row[f"near_identical_share_{number}"] = near.n_shared.sum() / max(genes.n_shared.sum(), 1)
        row[f"near_identical_kb_{number}"] = (near.end - near.start + 1).sum() / 1000
        # contigs on which nearly all genes are near-identical
        per_contig = genes.groupby("contig").near_identical.agg(["sum", "size"])
        wholly_shared = per_contig[per_contig["sum"] / per_contig["size"] >= WHOLLY_SHARED]
        row[f"wholly_shared_contigs_{number}"] = len(wholly_shared)
        row[f"wholly_shared_fraction_{number}"] = wholly_shared["sum"].sum() / max(len(near), 1)   # of the near-identical genes
        row[f"completeness_{number}"] = genomes.checkm_completeness[getattr(pair, f"genome_{number}")]
        row[f"contamination_{number}"] = genomes.checkm_contamination[getattr(pair, f"genome_{number}")]
    row["near_identical_share"] = max(row["near_identical_share_1"], row["near_identical_share_2"])
    row["near_identical_kb"] = max(row["near_identical_kb_1"], row["near_identical_kb_2"])
    row["wholly_shared_contigs"] = max(row["wholly_shared_contigs_1"], row["wholly_shared_contigs_2"])

    # the class of the pair (rules in the header)
    fraction_of_mags = [row[f"wholly_shared_fraction_{n}"] for n in (1, 2) if row[f"type_{n}"] != "isolate"]
    if row["near_identical_share"] < 0.5:
        row["class"] = "related genomes"
    elif fraction_of_mags:
        if max(fraction_of_mags) < 0.2:
            row["class"] = "MAG or SAG, unresolved"
        elif row["n_near_identical_orthologs"] >= 0.5 * row["n_orthologs"]:
            row["class"] = "MAG or SAG, close relatives with partial overlap"
        else:
            row["class"] = "MAG or SAG, foreign contigs"
    elif max(row["wholly_shared_fraction_1"], row["wholly_shared_fraction_2"]) >= 0.5:
        row["class"] = "isolates with wholly shared contigs"
    else:
        row["class"] = "isolates with embedded shared DNA"
    rows.append(row)

table = pd.DataFrame(rows).round(4)
table.to_csv("all_pairs/results/all_pairs_summary.tsv", sep="\t", index=False)

largest_difference = (table.fmh_dnds - table.fmh_dnds_in_paper).abs().max()
print(f"{len(table)} pairs; largest difference to the FracMinHash dN/dS of the manuscript's figure: {largest_difference:.1e}")
print(pd.crosstab(table["class"], table.fmh_dnds > 1).rename(columns={False: "at most 1", True: "above 1"}).to_string())
