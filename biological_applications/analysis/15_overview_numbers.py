#!/usr/bin/env python
"""Step 15. Numbers that the manuscript quotes about the 140 genome pairs as a whole, about the
effect of an identical fraction of sequence, and about the gene-removal experiment.

    usage:  python analysis/15_overview_numbers.py      (after steps 12 to 14 and after the pairs
                                                         listed in pairs/examined_pairs.tsv)

Inputs
    all_pairs/results/all_pairs_summary.tsv     one row per pair (step 12)
    all_pairs/results/gene_removal.tsv          the gene-removal experiment (step 13)
    all_pairs/results/translation_table_check.tsv   (step 14)
    pairs/examined_pairs.tsv                    the pairs examined in detail, and
    <their folder>/results/fmh_dnds.tsv, orthologs.tsv

The effect of an identical fraction
    FracMinHash dN/dS turns the nucleotide containment C_nt and the protein containment C_aa into
        dN/dS = (S/N) (1 - C_aa^(1/7)) / (C_aa^(1/7) - C_nt^(1/7)),    S/N = 0.77/2.23.
    If a fraction f of the coding sequence of two genomes is identical and the rest has diverged
    to the containments (c_nt, c_aa), then
        C_nt = f + (1 - f) c_nt        C_aa = f + (1 - f) c_aa.
    For each examined pair, c_nt and c_aa are the containments measured without the near-identical
    genes, and f is the value that reproduces the measured C_nt.

Outputs
    output/numbers/overview.tex                         LaTeX macros, prefix OV
    all_pairs/results/identical_fraction_model.tsv      per examined pair: f, predicted and observed dN/dS
"""
import numpy as np
import pandas as pd
from common import write_numbers, thousands, fmh_dnds_with_identical_fraction

NUMBER_WORDS = ["One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine"]
CLASS_MACROS = {      # class in the table -> macro name
    "related genomes": "ClassRelated",
    "MAG or SAG, foreign contigs": "ClassForeign",
    "MAG or SAG, close relatives with partial overlap": "ClassClose",
    "MAG or SAG, unresolved": "ClassMagInside",
    "isolates with wholly shared contigs": "ClassIsolateContigs",
    "isolates with embedded shared DNA": "ClassIsolateInside",
}

pairs = pd.read_csv("all_pairs/results/all_pairs_summary.tsv", sep="\t")
genomes = pd.read_csv("all_pairs/genomes.tsv", sep="\t")
above = pairs[pairs.fmh_dnds > 1]
below = pairs[pairs.fmh_dnds <= 1]
related = pairs[pairs["class"] == "related genomes"]
stayed = above[above.fmh_dnds_without_near_identical > 1]

# ------------------------------------------------------------------ all pairs
numbers = [
    ("Pairs", str(len(pairs)), "genome pairs in the hierarchical edge bundling figure"),
    ("Genomes", str(len(genomes)), "genomes in these pairs"),
    ("GenomesMagSag", str(int((genomes.genome_type != "isolate").sum())), "of those, metagenome-assembled or single-cell genomes"),
    ("PairsAboveOne", str(len(above)), "pairs with FracMinHash dN/dS above 1"),
    ("AboveOneNearDriven", str(int((above.near_identical_share >= 0.5).sum())),
     "of those, pairs in which at least half of the shared 21-mers lie in near-identical genes"),
    ("NearShareAbove", f"{100 * above.near_identical_share.median():.0f}", "median share of shared 21-mers in near-identical genes, pairs above 1 (percent)"),
    ("NearShareBelow", f"{100 * below.near_identical_share.median():.0f}", "the same for pairs at or below 1 (percent)"),
    ("MedianAbove", f"{above.fmh_dnds.median():.2f}", "median FracMinHash dN/dS of the pairs above 1"),
    ("MedianAboveWithout", f"{above.fmh_dnds_without_near_identical.median():.2f}", "the same pairs after removing near-identical genes"),
    ("FellBelowOne", str(int((above.fmh_dnds_without_near_identical <= 1).sum())), "pairs above 1 that fall to 1 or below after removal"),
    ("StayedAboveOne", str(len(stayed)), "pairs above 1 that stay above 1 after removal"),
    ("NoEstimateLeft", str(int(above.fmh_dnds_without_near_identical.isna().sum())), "pairs above 1 with no shared k-mers left after removal (no estimate)"),
    ("RelatedMedian", f"{related.fmh_dnds.median():.2f}", "median FracMinHash dN/dS of pairs of related genomes"),
    ("RelatedMedianWithout", f"{related.fmh_dnds_without_near_identical.median():.2f}", "the same after removing near-identical genes"),
    ("ConvMin", f"{pairs.conventional_median_dnds.min():.3f}", "smallest per-pair median of the alignment-based dN/dS"),
    ("ConvMax", f"{pairs.conventional_median_dnds.max():.2f}", "largest per-pair median of the alignment-based dN/dS"),
    ("Orthologs", thousands(pairs.n_orthologs.sum()), "orthologs aligned over all pairs"),
    ("MinContainment", f"{pairs.containment_nt.min():.3f}", "smallest nucleotide containment among the pairs of the figure"),
    ("MinAbove", f"{above.fmh_dnds.min():.2f}", "smallest FracMinHash dN/dS among the pairs above 1"),
    ("MaxBelow", f"{below.fmh_dnds.max():.2f}", "largest FracMinHash dN/dS among the pairs at or below 1"),
    # the same computation with every k-mer of the two genomes instead of a sketch (scale factor 1 instead of 500)
    ("EveryKmerAboveOne", str(int((above.fmh_dnds_every_kmer > 1).sum())), "pairs above 1 that are also above 1 when every k-mer is used"),
    ("EveryKmerBelowAboveOne", str(int((below.fmh_dnds_every_kmer > 1).sum())), "pairs at or below 1 that are above 1 when every k-mer is used"),
    ("EveryKmerDifferenceAbove", f"{(above.fmh_dnds - above.fmh_dnds_every_kmer).abs().median():.2f}", "median absolute difference between sketch and every-k-mer values, pairs above 1"),
    ("EveryKmerDifferenceBelow", f"{(below.fmh_dnds - below.fmh_dnds_every_kmer).abs().median():.2f}", "the same for pairs at or below 1"),
]
for name, macro in CLASS_MACROS.items():
    mine = pairs[pairs["class"] == name]
    numbers += [(macro, str(len(mine)), f"pairs in class '{name}'"),
                (macro + "AboveOne", str(int((mine.fmh_dnds > 1).sum())), "of those, with FracMinHash dN/dS above 1")]
for name, macro in CLASS_MACROS.items():
    mine = stayed[stayed["class"] == name]
    numbers.append(("Stayed" + macro, str(len(mine)), f"pairs that stay above 1 after removal, class '{name}'"))

# ------------------------------------------------------------------ MAG and SAG pairs with wholly shared contigs
foreign = pairs[pairs["class"] == "MAG or SAG, foreign contigs"]
close = pairs[pairs["class"] == "MAG or SAG, close relatives with partial overlap"]
ranks = foreign.lowest_shared_rank.value_counts()
numbers += [
    ("WhollyShared", str(len(foreign) + len(close)), "MAG/SAG pairs whose near-identical genes lie on wholly shared contigs"),
    ("WhollySharedAboveOne", str(int((foreign.fmh_dnds > 1).sum() + (close.fmh_dnds > 1).sum())), "of those, above 1"),
    ("ForeignSameGenus", str(int(ranks.get("genus", 0))), "foreign-contig pairs within one genus"),
    ("ForeignSameFamily", str(int(ranks.get("family", 0))), "foreign-contig pairs within one family (different genera)"),
    ("ForeignHigher", str(int(len(foreign) - ranks.get("genus", 0) - ranks.get("family", 0))), "foreign-contig pairs from different families or higher taxa"),
    ("ForeignCrossPhylum", str(int(ranks.get("domain", 0) + ranks.get("none", 0))), "of those, pairs from different phyla or different domains"),
    ("ForeignMedianContigs", f"{foreign.wholly_shared_contigs.median():.0f}", "foreign-contig pairs: median number of wholly shared contigs"),
    ("ForeignFewContigs", str(int((foreign.wholly_shared_contigs <= 2).sum())), "foreign-contig pairs with only one or two wholly shared contigs"),
    ("ForeignMedianKb", f"{foreign.near_identical_kb.median():.0f}", "foreign-contig pairs: median length of near-identical genes (kb)"),
    ("ForeignNearOrthologs", f"{100 * (foreign.n_near_identical_orthologs / foreign.n_orthologs).median():.0f}", "foreign-contig pairs: median percent of orthologs that are near-identical"),
    ("ForeignBackgroundIdentity", f"{foreign.background_aa_identity.median():.0f}", "foreign-contig pairs: median amino-acid identity of the other orthologs (percent)"),
    ("CloseNearOrthologs", f"{100 * (close.n_near_identical_orthologs / close.n_orthologs).median():.0f}", "close-relative pairs: median percent of orthologs that are near-identical"),
    ("CloseSameGenus", str(int((close.lowest_shared_rank == "genus").sum())), "close-relative pairs within one genus"),
    ("CloseCompleteness", f"{close[['completeness_1', 'completeness_2']].min(axis=1).median():.0f}", "close-relative pairs: median CheckM completeness of the less complete genome (percent)"),
]
# the MAGs and SAGs that carry the wholly shared contigs of the foreign-contig pairs, and their CheckM contamination
carriers = {}
for row in foreign.itertuples():
    for accession, kind, fraction, contamination in ((row.genome_1, row.type_1, row.wholly_shared_fraction_1, row.contamination_1),
                                                    (row.genome_2, row.type_2, row.wholly_shared_fraction_2, row.contamination_2)):
        if kind != "isolate" and fraction >= 0.2:
            carriers[accession] = contamination
contamination = pd.Series(carriers)
# share of a MAG's or SAG's near-identical genes that lies on fully shared contigs (the larger value of a pair)
on_contigs = [max(fraction for fraction, kind in ((row.wholly_shared_fraction_1, row.type_1), (row.wholly_shared_fraction_2, row.type_2)) if kind != "isolate")
              for row in foreign.itertuples()]
numbers.append(("ForeignShareOnContigs", f"{100 * pd.Series(on_contigs).median():.0f}",
                "foreign-contig pairs: median percent of the near-identical genes of the MAG or SAG that lie on fully shared contigs"))
contamination.rename("checkm_contamination").rename_axis("accession").to_csv("all_pairs/results/carrier_genomes.tsv", sep="\t")
numbers += [
    ("CarrierGenomes", str(len(contamination)), "MAGs and SAGs that carry the wholly shared contigs of foreign-contig pairs"),
    ("CarrierContamination", f"{contamination.median():.2f}", "their median CheckM contamination (percent)"),
    ("CarrierPass", str(int((contamination <= 5).sum())), "of those genomes, with CheckM contamination of at most 5 percent"),
]
isolate_contigs = pairs[(pairs["class"] == "isolates with wholly shared contigs") & (pairs.wholly_shared_contigs >= 10)]
numbers += [
    ("IsolateContaminated", str(len(isolate_contigs)), "isolate pairs with at least ten wholly shared contigs"),
    ("IsolateContigsMin", str(int(isolate_contigs.wholly_shared_contigs.min())), "fewest wholly shared contigs among them"),
    ("IsolateContigsMax", str(int(isolate_contigs.wholly_shared_contigs.max())), "most wholly shared contigs among them"),
]

# ------------------------------------------------------------------ the effect of an identical fraction, for the examined pairs
examined = pd.read_csv("pairs/examined_pairs.tsv", sep="\t")
wide = np.linspace(0, 0.5, 5001)
model_rows = []
for pair in examined.itertuples():
    fmh = pd.read_csv(f"{pair.folder}/results/fmh_dnds.tsv", sep="\t").set_index("genes")
    orthologs = pd.read_csv(f"{pair.folder}/results/orthologs.tsv", sep="\t")
    c_nt, c_aa = fmh.loc["without_near_identical", ["containment_nt", "containment_aa"]]
    row = dict(number=pair.number, name=pair.name, observed_fmh_dnds=fmh.loc["all", "fmh_dnds"],
               background_fmh_dnds=fmh.loc["without_near_identical", "fmh_dnds"], background_containment_nt=c_nt, background_containment_aa=c_aa,
               background_aa_identity=orthologs[~orthologs.near_identical].aa_identity.median())
    if c_nt > 0:        # a pair with no shared k-mers outside the near-identical genes has no background to model
        row["identical_fraction"] = (fmh.loc["all", "containment_nt"] - c_nt) / (1 - c_nt)
        row["predicted_fmh_dnds"] = fmh_dnds_with_identical_fraction(row["identical_fraction"], c_nt, c_aa)
        curve = fmh_dnds_with_identical_fraction(wide, c_nt, c_aa)
        row["fraction_at_which_estimate_exceeds_1"] = wide[np.argmax(curve > 1)] if (curve > 1).any() else float("nan")
    model_rows.append(row)
model = pd.DataFrame(model_rows)
model.to_csv("all_pairs/results/identical_fraction_model.tsv", sep="\t", index=False)
for row in model.dropna(subset=["identical_fraction"]).itertuples():
    word = NUMBER_WORDS[row.number - 1]
    crossing = "not reached" if np.isnan(row.fraction_at_which_estimate_exceeds_1) else f"{100 * row.fraction_at_which_estimate_exceeds_1:.1f}"
    numbers += [(f"Background{word}", f"{row.background_fmh_dnds:.2f}", f"pair {row.number} ({row.name}): FracMinHash dN/dS without near-identical genes"),
                (f"Observed{word}", f"{row.observed_fmh_dnds:.2f}", f"pair {row.number}: observed FracMinHash dN/dS"),
                (f"Predicted{word}", f"{row.predicted_fmh_dnds:.2f}", f"pair {row.number}: value predicted for the observed identical fraction"),
                (f"IdenticalPercent{word}", f"{100 * row.identical_fraction:.1f}", f"pair {row.number}: identical fraction that reproduces the nucleotide containment (percent)"),
                (f"CrossingPercent{word}", crossing, f"pair {row.number}: identical fraction at which the estimate exceeds 1 (percent)")]
modelled = model.dropna(subset=["identical_fraction"])
high, low = modelled[modelled.observed_fmh_dnds > 1], modelled[modelled.observed_fmh_dnds <= 1]
numbers += [("ModelPairsAbove", str(len(high)), "examined pairs above 1 to which the model was applied"),
            ("ModelPairsBelow", str(len(low)), "examined pairs below 1 to which the model was applied"),
            ("ModelLargestDeviation", f"{(modelled.predicted_fmh_dnds - modelled.observed_fmh_dnds).abs().max():.2f}", "largest absolute difference between predicted and observed values"),
            ("CrossingAboveMin", f"{100 * high.fraction_at_which_estimate_exceeds_1.min():.1f}", "pairs above 1: smallest identical fraction at which the estimate exceeds 1 (percent)"),
            ("CrossingAboveMax", f"{100 * high.fraction_at_which_estimate_exceeds_1.max():.1f}", "pairs above 1: largest such fraction (percent)"),
            ("CrossingBelowMin", f"{100 * low.fraction_at_which_estimate_exceeds_1.min():.1f}", "pairs below 1: smallest identical fraction at which the estimate would exceed 1 (percent)"),
            ("CrossingBelowMax", f"{100 * low.fraction_at_which_estimate_exceeds_1.max():.1f}", "pairs below 1: largest such fraction (percent)"),
            ("IdentityAboveMin", f"{high.background_aa_identity.min():.0f}", "pairs above 1: lowest median amino-acid identity of the diverged orthologs (percent)"),
            ("IdentityAboveMax", f"{high.background_aa_identity.max():.0f}", "pairs above 1: highest median amino-acid identity of the diverged orthologs (percent)"),
            ("IdentityBelowMin", f"{low.background_aa_identity.min():.0f}", "pairs below 1: lowest median amino-acid identity of the diverged orthologs (percent)"),
            ("IdentityBelowMax", f"{low.background_aa_identity.max():.0f}", "pairs below 1: highest median amino-acid identity of the diverged orthologs (percent)")]

# ------------------------------------------------------------------ the gene-removal experiment
removal = pd.read_csv("all_pairs/results/gene_removal.tsv", sep="\t")
for (pair, mine), word in zip(removal.groupby("pair", sort=False), ["One", "Two"]):
    means = mine.groupby("kept_percent").fmh_dnds.mean()
    for kept, kept_word in {100: "All", 90: "Ninety", 70: "Seventy", 50: "Half"}.items():
        numbers.append((f"Removal{word}Keep{kept_word}", f"{means[kept]:.2f}",
                        f"gene-removal experiment, {pair}: mean FracMinHash dN/dS when each genome keeps {kept} percent of its genes"))

# ------------------------------------------------------------------ the pair whose genes were called with translation table 4
check = pd.read_csv("all_pairs/results/translation_table_check.tsv", sep="\t").set_index("translation_table")
numbers += [("TableFourFmhDnds", f"{check.fmh_dnds[4]:.2f}", "Neisseria pair: FracMinHash dN/dS with the GTDB gene calls of GCA_923184405.1 (translation table 4)"),
            ("TableElevenFmhDnds", f"{check.fmh_dnds[11]:.2f}", "Neisseria pair: FracMinHash dN/dS when the genes of GCA_923184405.1 are called with table 11"),
            ("TableFourGenes", thousands(check.genes[4]), "genes of GCA_923184405.1 predicted with table 4"),
            ("TableElevenGenes", thousands(check.genes[11]), "genes of GCA_923184405.1 predicted with table 11")]

write_numbers("output/numbers/overview.tex", "OV", numbers)
