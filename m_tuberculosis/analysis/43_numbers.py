"""Every number the vignette quotes, as LaTeX macros, from the tool-based results.

    usage:  python analysis/43_numbers.py

Each macro is a \\providecommand with a comment saying what it is. Macro names cannot contain
digits, so scale factors are spelled out (Two, Three, Five, Ten, Hundred, Thousand).
out: results/numbers/mtb.tex
"""
import json
import os
import re

import numpy as np
import pandas as pd

MTB = os.getcwd()
OUT = f"{MTB}/results/numbers"
os.makedirs(OUT, exist_ok=True)
SPELLED = {1: "One", 2: "Two", 3: "Three", 5: "Five", 10: "Ten", 100: "Hundred", 1000: "Thousand"}
GROUP_NAMES = {"all_cds": "Whole", "essential": "Essential", "nonessential": "Nonessential",
               "pe_ppe": "Peppe", "pe_ppe_hq": "Peppehq"}
macros = []


def add(name, value, comment):
    macros.append((name, value, comment))


def number(value, decimals=0):
    if decimals == 0:
        return f"{round(float(value)):,}".replace(",", "{,}")
    text = f"{float(value):.{decimals}f}"
    return text[1:] if text.startswith("-") and float(text) == 0 else text


def read(path, **kwargs):
    return pd.read_csv(path, sep="\t", **kwargs) if os.path.exists(path) else None


# ---- the collection
collection = read(f"{MTB}/results/data_overlap/collection.tsv").iloc[0]
inputs = read(f"{MTB}/results/inputs/isolates.tsv")
isolates = inputs[inputs.source != "reference"]
add("mtbRunsGan", number(collection.runs), "distinct runs in Gan et al. Table S1")
add("mtbWithAssembly", number(collection.runs_with_assembly), "runs with an AllTheBacteria 0.2 assembly")
add("mtbAssemblyIsolates", number((isolates.source == "assembly").sum()), "isolates from assemblies")
add("mtbReadIsolates", number((isolates.source == "reads").sum()), "isolates from reads only (no assembly)")
add("mtbAllIsolates", number(len(isolates)), "isolates in total")
add("mtbMeanSubstitutions", number(isolates.substitutions.mean()), "mean substitutions per isolate")
coverage = read(f"{MTB}/results/inputs/coverage.tsv")
add("mtbMeanCoverage", number(100 * coverage.covered_bases.mean() / 4411532, 1),
    "mean percentage of the reference covered by an assembly")
groups_table = read(f"{MTB}/results/groups/summary.tsv").set_index("group")
for group, name in GROUP_NAMES.items():
    if group in groups_table.index:
        add(f"mtbGenes{name}", number(groups_table.loc[group, "genes"]), f"coding sequences, {group}")
        add(f"mtbBases{name}", number(groups_table.loc[group, "coding_bases"]), f"coding bases, {group}")

# ---- draft assemblies (the tool's main mode, step 12)
smoke = read(f"{MTB}/results/smoke/summary.tsv")
exact = smoke[smoke.scale == 1].iloc[0]
add("mtbSmokeIsolates", number(len(read(f"{MTB}/results/smoke/isolates.tsv"))), "isolates, draft-assembly test")
add("mtbSmokePairs", number(exact.pairs), "pairs, draft-assembly test")
add("mtbAssemblyMedianDnds", number(exact.median_dnds, 2), "median dN/dS from draft assemblies, every k-mer")
add("mtbAssemblyUnusable", number(100 * (1 - exact.fraction_usable)), "percentage of those pairs with no usable value")
add("mtbAssemblyCnt", number(1 - exact.median_C_nt, 4), "1 - C_nt, draft assemblies")
add("mtbAssemblyCaa", number(1 - exact.median_C_aa, 4), "1 - C_aa, draft assemblies")

# ---- the tool's runs
trends = f"{MTB}/results/trends_tool"
groups = read(f"{trends}/groups.tsv").set_index("group")
for group, name in GROUP_NAMES.items():
    if group in groups.index:
        row = groups.loc[group]
        add(f"mtbMedian{name}", number(row["median"], 3), f"median dN/dS, {group}, scale factor 1, every pair")
        add(f"mtbQuartiles{name}", f"{number(row.q1, 3)}--{number(row.q3, 3)}", f"quartiles, {group}")
        add(f"mtbAboveOne{name}", number(100 * row.above_one, 1), f"percentage of pairs above 1, {group}")
whole = groups.loc["all_cds"]
add("mtbPairs", number(whole.pairs), "pairs of assembly-derived isolates")
add("mtbPairsMillion", number(whole.pairs / 1e6), "the same, in millions")
add("mtbBelowOne", number(100 * whole.below_one, 2), "percentage of pairs below 1, whole coding genome")
categories = read(f"{trends}/whole_genome.tsv").set_index("category")
for category, name in (("same sublineage", "SameSublineage"),
                       ("same lineage, different sublineage", "SameLineage"),
                       ("different lineage", "DifferentLineage")):
    add(f"mtbMedian{name}", number(categories.loc[category, "median"], 3), f"median dN/dS, {category}")
    add(f"mtbPairs{name}", number(categories.loc[category, "pairs"]), f"pairs, {category}")

# what the figures show, written by 42_figures.py
figure_parameters = read(f"{MTB}/results/figures/parameters.tsv").set_index("name").value.astype(float)

paired = read(f"{trends}/paired.tsv")
for _, row in paired.iterrows():
    a, b = GROUP_NAMES[row.group_a], GROUP_NAMES[row.group_b]
    add(f"mtbShare{b}Over{a}", number(100 * row.share_b_greater, 1), f"percentage of pairs where {row.group_b} > {row.group_a}")
    add(f"mtbShift{b}Over{a}", number(row.median_difference, 3), f"median of {row.group_b} minus {row.group_a}")
add("mtbPairedPairsMillion", number(paired.pairs.min() / 1e6, 1), "pairs in the shared subsample, millions")
tested = paired[((paired.group_a == "essential") & (paired.group_b == "nonessential")) |
                ((paired.group_a == "nonessential") & (paired.group_b == "pe_ppe_hq"))]
add("mtbWilcoxonPairs", number(int(tested.wilcoxon_pairs.max())), "pairs in each Wilcoxon signed-rank test")
# a p-value of 0.0 is below the smallest positive double (about 1e-308), so 1e-300 bounds it
largest_p = tested.wilcoxon_p.max()
add("mtbWilcoxonPExponent", number(int(np.floor(np.log10(largest_p))) + 1 if largest_p > 0 else -300),
    "p < 10^this for both Wilcoxon tests of the ordering")

scales = read(f"{trends}/scales.tsv")
for _, row in scales.iterrows():
    if row.scale == 1:
        continue
    name = f"{GROUP_NAMES[row.group]}Scale{SPELLED[int(row.scale)]}"
    add(f"mtbSpearman{name}", number(row.spearman_with_scale_1, 2), f"rank correlation with scale 1, {row.group}, scale {row.scale}")
    add(f"mtbMedian{name}", number(row["median"], 3), f"median dN/dS, {row.group}, scale {row.scale}")
    add(f"mtbAboveOne{name}", number(100 * row.above_one, 1), f"percentage above 1, {row.group}, scale {row.scale}")
ordering = read(f"{trends}/ordering.tsv")
for _, row in ordering.iterrows():
    name = f"{GROUP_NAMES[row.lower]}Below{GROUP_NAMES[row.higher]}Scale{SPELLED[int(row.scale)]}"
    add(f"mtbHolds{name}", number(100 * row.share_holding, 1),
        f"percentage of pairs where {row.lower} < {row.higher}, scale {row.scale}")

tail = read(f"{trends}/divergence_tail_summary.tsv")
if tail is not None:
    tail = tail.iloc[0]
    add("mtbTailThreshold", number(100 * tail.threshold, 2), "per-nucleotide divergence above which pairs are rare (percent)")
    add("mtbTailShare", number(100 * tail.tail_share, 3), "percentage of pairs above that divergence")
    add("mtbTailIsolates", number(tail.outlier_isolates), "isolates that are in that tail for most of their pairs")
    add("mtbTailSubstitutions", number(tail.outlier_median_substitutions), "their median number of substitutions")
    add("mtbTailWithOutlier", number(100 * tail.tail_share_with_outlier), "percentage of tail pairs involving one of them")
divergence = read(f"{trends}/divergence.tsv")
if divergence is not None:
    shown = divergence[divergence.pairs >= figure_parameters["divergence_min_pairs"]]
    add("mtbDivergenceLowMedian", number(shown.iloc[0]["median"], 2), "median dN/dS of the least divergent bin shown")

every = read(f"{trends}/every_isolate.tsv")
if every is not None:
    every = every.set_index("category").loc["all pairs"]
    run = json.load(open(f"{MTB}/results/tool/all_cds_every_isolate/scaled1/run.json"))
    add("mtbFullPairs", number(every.pairs), "pairs over every isolate")
    add("mtbFullPairsBillion", number(every.pairs / 1e9, 2), "the same, in billions")
    add("mtbFullMedian", number(every["median"], 3), "median dN/dS over every isolate")
    add("mtbFullPairsMinutes", number(run["pairs_seconds"] / 60, 1), "minutes for the pairs, every isolate")
    add("mtbFullLoadMinutes", number(run["load_seconds"] / 60, 1), "minutes to read the sketches, every isolate")
    add("mtbCores", number(run["cores"]), "cores used")

timings = read(f"{trends}/timings.tsv")
row = timings[(timings.run == "all_cds") & (timings.scale == 1)].iloc[0]
add("mtbPairsMinutes", number(row.pairs_seconds / 60, 1), "minutes for all pairs, whole coding genome, scale 1")
add("mtbLoadMinutes", number(row.load_seconds / 60, 1), "minutes to read the sketches and form the differences")
add("mtbCoreSecondsPerMillion", number(row.core_seconds_per_million_pairs, 1),
    "core-seconds per million pairs, whole coding genome, scale 1")

full_run = f"{MTB}/results/tool/all_cds_every_isolate/scaled1/run.json"
sketching = read(f"{MTB}/results/tool/sketching.tsv")
if os.path.exists(full_run) and sketching is not None:
    report = json.load(open(full_run))
    sketch = sketching.set_index("run").loc["all_cds_every_isolate"]
    add("mtbSketchMinutes", number(sketch.wall_seconds / 60), "minutes the tool took to sketch every genome, scale factor 1")
    add("mtbSketchGenomes", number(sketch.genomes), "genomes sketched (the isolates and the reference)")
    add("mtbSketchShards", number(sketch.shards), "shards sketched side by side, one core each")
    pairs_core_hours = (report["load_seconds"] + report["pairs_seconds"]) * report["cores"] / 3600
    add("mtbToolCoreHours", number(sketch.core_hours + pairs_core_hours),
        "core-hours from FASTA files to every pair's dN/dS, every isolate (wall time times cores)")
    add("mtbTotalMinutes", number((sketch.wall_seconds + report["load_seconds"] + report["pairs_seconds"]) / 60),
        "minutes from FASTA files to every pair's dN/dS, every isolate")

# ---- seqinr on the same sequences against the tool's main mode
conventional = read(f"{MTB}/results/conventional_tool/summary.tsv").set_index("group")
for group, name in GROUP_NAMES.items():
    if group in conventional.index:
        row = conventional.loc[group]
        add(f"mtbSeqinrMedian{name}", number(row.median_seqinr, 3), f"seqinr median, {group}")
        add(f"mtbToolMedian{name}", number(row.median_fmh, 3), f"the tool's median on those pairs, {group}")
        add(f"mtbSpearmanSeqinr{name}", number(row.spearman, 2), f"rank correlation, tool and seqinr, {group}")
row = conventional.loc["all_cds"]
add("mtbSeqinrIsolates", number(row.isolates), "isolates compared with seqinr")
add("mtbSeqinrPairs", number(row.pairs), "pairs compared with seqinr")
add("mtbSeqinrSeconds", number(row.seqinr_seconds, 1), "seconds seqinr took for those pairs, one core")
per_pair = row.seqinr_seconds / row.pairs
add("mtbSeqinrPairsPerSecond", number(1 / per_pair), "pairs per second, seqinr, one core, whole coding genome")
n_all = len(isolates)
add("mtbSeqinrAllCoreYears", number(per_pair * n_all * (n_all - 1) / 2 / (3600 * 24 * 365), 1),
    "core-years seqinr would need for every pair of the collection")
add("mtbSeqinrFiveThousandCoreDays", number(per_pair * 5000 * 4999 / 2 / (3600 * 24), 1),
    "core-days seqinr would need for every pair of 5,000 isolates")

# ---- the constant, the differing hashes, and reads against assemblies
site_ratio = read(f"{MTB}/results/site_ratio.tsv").set_index("group")
add("mtbSiteRatioWhole", number(site_ratio.loc["all_cds", "constant_for_these_sequences"], 3),
    "synonymous over non-synonymous sites in these coding sequences")
add("mtbConstant", number(site_ratio.loc["all_cds", "constant_used"], 3), "the constant the estimator uses, 0.77/2.23")
for _, row in scales[scales.scale == 1].iterrows():
    add(f"mtbDifferingProtein{GROUP_NAMES[row.group]}", number(row.differing_protein_hashes),
        f"median protein hashes that differ between the two isolates of a pair, {row.group}, scale factor 1")
    add(f"mtbDifferingDna{GROUP_NAMES[row.group]}", number(row.differing_dna_hashes),
        f"median nucleotide hashes that differ between the two isolates of a pair, {row.group}, scale factor 1")
for _, row in scales[scales.scale > 1].iterrows():
    add(f"mtbDifferingProtein{GROUP_NAMES[row.group]}Scale{SPELLED[int(row.scale)]}", number(row.differing_protein_hashes),
        f"median differing protein hashes per pair, {row.group}, scale factor {row.scale}")
reads = f"{MTB}/results/reads_vs_assembly_tool/summary.txt"
if os.path.exists(reads):
    text = open(reads).read()
    add("mtbReadsAssemblyIsolates", re.search(r"isolates with both kinds of input: (\d+)", text).group(1),
        "isolates with both an assembly and reads")
    add("mtbReadsAssemblySpearman", number(float(re.search(r"Spearman correlation: ([0-9.]+)", text).group(1)), 2),
        "rank correlation of dN/dS from reads and from assemblies, same pairs")
    add("mtbReadsAssemblyShift", number(float(re.search(r"reads minus assemblies: (-?[0-9.]+)", text).group(1)), 3),
        "median difference of dN/dS, reads minus assemblies")

# ---- how many whole coding genomes main mode's single sketch file holds before it passes 4 GB
sketch_file = read(f"{MTB}/results/checks/main_mode_sketch.tsv").iloc[0]
add("mtbMainModeLimitGenomes", number(int(2 ** 32 // sketch_file.bytes_per_genome)),
    "whole coding genomes at scale factor 1 that fit in main mode's sketch file below 4 GiB (2^32 bytes)")
add("mtbSketchMegabytesPerGenome", number(sketch_file.bytes_per_genome / 1e6), "megabytes of sketch per whole coding genome, scale factor 1")

# ---- what the figures show (written by 42_figures.py)
add("mtbSeqinrPanelLimit", number(figure_parameters["seqinr_limit"], 1), "upper limit of both axes of the seqinr panel")
add("mtbSeqinrPanelOutside", number(int(figure_parameters["seqinr_pairs_outside"])), "whole-genome pairs beyond the seqinr panel's axes")
add("mtbDivergenceMinPairs", number(int(figure_parameters["divergence_min_pairs"])), "fewest pairs in a divergence bin that is drawn")
add("mtbHistogramBinWidth", number(figure_parameters["histogram_bin_width"], 3), "width of the bins of the distribution panel")
add("mtbAnchoredPoints", number(int(figure_parameters["anchored_points"])), "anchored pairs drawn in the containment panel")

# ---- checks: main mode against reference mode (step 31, summarised by step 31b)
check = read(f"{MTB}/results/checks/reference_mode.tsv")
assert check.identical.all(), "the check found differences"
add("mtbEquivalenceComparisons", number(len(check)), "comparisons of main and reference mode (group x scale factor), all identical")
add("mtbEquivalencePairs", number(check.pairs.sum()), "pairs compared over those comparisons")
add("mtbEquivalencePairsSegments", number(check[check.group == "pe_ppe_hq"].pairs.max()), "pairs per comparison, PE/PPE segments (1,000 isolates)")
add("mtbEquivalencePairsWhole", number(check[check.group == "all_cds"].pairs.max()), "pairs per comparison, whole coding genome (80 isolates)")
prototype = read(f"{trends}/prototype.tsv")
if prototype is not None:
    add("mtbPrototypeDifferentBins", number(prototype.pairs_in_different_bins.sum()),
        "pairs whose value falls in a different 0.005-wide bin in the earlier prototype")

lines = ["% Numbers of the M. tuberculosis vignette, written by analysis/43_numbers.py from the",
         "% FracMinHash dN/dS tool's output. Do not edit by hand: rerun the script instead.", ""]
for name, value, comment in macros:
    lines += [f"% {comment}", f"\\providecommand{{\\{name}}}{{{value}}}"]
open(f"{OUT}/mtb.tex", "w").write("\n".join(lines) + "\n")
bad = [name for name, _, _ in macros if re.search(r"[^A-Za-z]", name)]
assert not bad, f"macro names with characters LaTeX does not allow: {bad}"
print(f"wrote {len(macros)} macros to {OUT}/mtb.tex")
