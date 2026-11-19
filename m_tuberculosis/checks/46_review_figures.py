"""Adversarial review of the vignette's figures: every plotted quantity recomputed from the tool's raw
per-pair output, without the code of steps 33 and 40-42, and compared with what the panels draw and
the captions and text state.

    usage:  python checks/46_review_figures.py [--jobs N]

Independent of the pipeline in three ways:
  * quantiles are found by a different algorithm (a fine histogram locates the bin of each order
    statistic, a second pass collects the values in those bins) from the raw parquet parts;
  * pair-by-pair comparisons use a different subsample of pairs (selected by a CRC32 hash of the two
    isolate names, not by their positions), so agreement is up to sampling error only;
  * seqinr's kaks is run again on two-sequence alignments of randomly chosen pairs.

out: results/review/checks.tsv   one row per check: figure panel, what is checked, the value drawn or
                                  stated, the value recomputed, the difference, the tolerance, pass/fail
     results/review/seqinr_rerun.tsv   the re-run seqinr pairs
"""
import argparse
import glob
import json
import os
import random
import re
import subprocess
import zlib
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import stats

MTB = os.getcwd()
OUT = f"{MTB}/results/review"
REFERENCE = "ancestral_reference"
GROUPS = ["essential", "all_cds", "nonessential", "pe_ppe", "pe_ppe_hq"]
SCALES = [1, 2, 3, 5, 10, 100, 1000]
QUANTILES = {"p5": 0.05, "q1": 0.25, "median": 0.5, "q3": 0.75, "p95": 0.95}
EDGES = np.concatenate([[-np.inf], np.linspace(-10, 20, 150_001), [np.inf]])   # bins 0.0002 wide
SAMPLE_MODULUS, SAMPLE_REMAINDER = 64, 17
EXACT, SAMPLED = 1e-9, 0.01      # tolerances: recomputed exactly, or on a different subsample
checks = []


def check(panel, what, stated, recomputed, tolerance):
    difference = abs(float(stated) - float(recomputed))
    checks.append(dict(panel=panel, check=what, stated=stated, recomputed=recomputed, difference=difference,
                       tolerance=tolerance, verdict="pass" if difference <= tolerance else "FAIL"))
    print(f"  {checks[-1]['verdict']:4}  {panel:8} {what}: stated {stated}, recomputed {recomputed}", flush=True)


def initialise(state):
    pa.set_cpu_count(1)              # one thread per worker process
    pa.set_io_thread_count(1)
    globals().update(state)


def codes(column):
    """Global index of each name in a dictionary-encoded column; -1 for the reference."""
    parts = []
    for chunk in column.chunks:
        mapping = np.array([INDEX.get(n, -1) for n in chunk.dictionary.to_pylist()], dtype=np.int64)
        parts.append(mapping[chunk.indices.to_numpy()])
    return np.concatenate(parts)


def read_part(path):
    table = pq.read_table(path, columns=["query_name", "match_name", "DNA_max_Cfrac", "AA_max_Cfrac", "dN/dS"],
                          use_threads=False)
    a, b = codes(table.column("query_name")), codes(table.column("match_name"))
    keep = (a >= 0) & (b >= 0)
    low, high = np.minimum(a, b)[keep], np.maximum(a, b)[keep]
    return (low * len(INDEX) + high, table.column("DNA_max_Cfrac").to_numpy()[keep],
            table.column("AA_max_Cfrac").to_numpy()[keep], table.column("dN/dS").to_numpy()[keep])


def first_pass(path):
    """Counts, a fine histogram of the finite values, this review's own subsample, and the rows of
    any pairs asked for by name."""
    key, c_nt, c_aa, dnds = read_part(path)
    finite = np.isfinite(dnds)
    low, high = key // len(INDEX), key % len(INDEX)
    sampled = (NAME_HASH[low] + NAME_HASH[high]) % SAMPLE_MODULUS == SAMPLE_REMAINDER
    wanted = np.isin(key, WANTED) if len(WANTED) else np.zeros(len(key), bool)
    return dict(pairs=len(key), finite=int(finite.sum()), below_one=int((dnds[finite] < 1).sum()),
                histogram=np.histogram(dnds[finite], bins=EDGES)[0],
                sample=(key[sampled], dnds[sampled]),
                wanted=pd.DataFrame(dict(key=key[wanted], C_nt=c_nt[wanted], C_aa=c_aa[wanted], dnds=dnds[wanted])))


def second_pass(path):
    """The finite values that fall in the bins holding the order statistics wanted, and, if asked,
    the values of the pairs in chosen divergence bins."""
    _, c_nt, _, dnds = read_part(path)
    finite = np.isfinite(dnds)
    values, c_nt = dnds[finite], c_nt[finite]
    in_bins = np.isin(np.searchsorted(EDGES, values, side="right") - 1, TARGET_BINS)
    by_divergence = {}
    if len(DIVERGENCE_BINS):
        divergence = 1.0 - np.power(np.clip(c_nt, 0, None), 1.0 / 21)
        for low, high in DIVERGENCE_BINS:
            by_divergence[(low, high)] = values[(divergence >= low) & (divergence < high)]
    return values[in_bins], by_divergence


def order_statistics(histogram, collected, ranks):
    """The values of the given ranks (0-based) among all finite values, from the histogram and the
    values collected from the bins holding them."""
    cumulative = np.cumsum(histogram)
    result = {}
    for rank in ranks:
        b = int(np.searchsorted(cumulative, rank, side="right"))
        below = cumulative[b - 1] if b else 0
        result[rank] = np.sort(collected[b])[rank - below]
    return result


def raw_run(run, scale, jobs, names, wanted=(), divergence_bins=()):
    """Everything this review needs from one run of the tool, in two passes over its parquet parts."""
    folder = f"{MTB}/work/tool/{run}/reference_scaled{scale}"
    parts = sorted(glob.glob(f"{folder}/pairs/part-*.parquet"))
    state = dict(INDEX={n: i for i, n in enumerate(names)},
                 NAME_HASH=np.array([zlib.crc32(n.encode()) for n in names], dtype=np.int64),
                 WANTED=np.array(sorted(wanted), dtype=np.int64), TARGET_BINS=np.array([], dtype=np.int64),
                 DIVERGENCE_BINS=())
    totals, histogram, sample_keys, sample_values, wanted_rows = dict(pairs=0, finite=0, below_one=0), 0, [], [], []
    with Pool(jobs, initializer=initialise, initargs=(state,)) as pool:
        for result in pool.imap_unordered(first_pass, parts, chunksize=1):
            for name in totals:
                totals[name] += result[name]
            histogram = histogram + result["histogram"]
            sample_keys.append(result["sample"][0])
            sample_values.append(result["sample"][1])
            wanted_rows.append(result["wanted"])
    n = totals["finite"]
    ranks = {}
    for name, q in QUANTILES.items():               # numpy's default: linear between two order statistics
        position = (n - 1) * q
        ranks[name] = (int(np.floor(position)), min(int(np.floor(position)) + 1, n - 1), position - np.floor(position))
    needed = sorted({r for low, high, _ in ranks.values() for r in (low, high)})
    cumulative = np.cumsum(histogram)
    bins = sorted({int(np.searchsorted(cumulative, r, side="right")) for r in needed})
    assert 0 not in bins and len(EDGES) - 2 not in bins, "an order statistic fell in an open-ended bin"
    state.update(TARGET_BINS=np.array(bins, dtype=np.int64), DIVERGENCE_BINS=tuple(divergence_bins))
    collected = {b: [] for b in bins}
    by_divergence = {d: [] for d in divergence_bins}
    with Pool(jobs, initializer=initialise, initargs=(state,)) as pool:
        for values, divergence_values in pool.imap_unordered(second_pass, parts, chunksize=1):
            which = np.searchsorted(EDGES, values, side="right") - 1
            for b in bins:
                collected[b].append(values[which == b])
            for d, v in divergence_values.items():
                by_divergence[d].append(v)
    collected = {b: np.concatenate(v) for b, v in collected.items()}
    found = order_statistics(histogram, collected, needed)
    quantiles = {name: found[low] + fraction * (found[high] - found[low]) for name, (low, high, fraction) in ranks.items()}
    sample = pd.Series(np.concatenate(sample_values), index=np.concatenate(sample_keys))
    return dict(totals=totals, quantiles=quantiles, sample=sample, wanted=pd.concat(wanted_rows),
                by_divergence={d: np.concatenate(v) for d, v in by_divergence.items()}, folder=folder)


def macros():
    text = open(f"{MTB}/results/numbers/mtb.tex").read()
    return dict(re.findall(r"\\providecommand\{\\(\w+)\}\{((?:[^{}]|\{[^{}]*\})*)\}", text))


def as_number(macro):
    return float(macro.replace("{,}", "").replace("$", ""))


def rerun_seqinr(task):
    """seqinr's kaks on one pair, as a two-sequence alignment."""
    group, a, b, folder = task
    path = f"{folder}/{group}_{a}_{b}.fasta"
    with open(path, "w") as out:
        for name in (a, b):
            sequence = "".join(line.strip() for line in open(f"{MTB}/work/fasta/{group}/{name}.ffn") if not line.startswith(">"))
            out.write(f">{name}\n{sequence}\n")
    script = ("suppressMessages(library(seqinr)); r <- kaks(read.alignment(commandArgs(TRUE)[1], format='fasta'));"
              " cat(as.numeric(r$ka), as.numeric(r$ks))")
    output = subprocess.run(["Rscript", "-e", script, path], capture_output=True, text=True, check=True).stdout
    os.remove(path)
    ka, ks = (float(x) for x in output.split())
    return dict(group=group, A=a, B=b, ka=ka, ks=ks)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--jobs", type=int, default=200)
    args = parser.parse_args()
    os.makedirs(OUT, exist_ok=True)
    stated = macros()
    trends = f"{MTB}/results/trends_tool"
    groups = pd.read_csv(f"{trends}/groups.tsv", sep="\t").set_index("group")
    scales = pd.read_csv(f"{trends}/scales.tsv", sep="\t").set_index(["group", "scale"])
    ordering = pd.read_csv(f"{trends}/ordering.tsv", sep="\t")
    divergence_drawn = pd.read_csv(f"{trends}/divergence.tsv", sep="\t")
    parameters = pd.read_csv(f"{MTB}/results/figures/parameters.tsv", sep="\t").set_index("name").value
    names = sorted(line.strip() for line in open(f"{MTB}/work/fasta/assembly_names.txt") if line.strip())
    index = {n: i for i, n in enumerate(names)}

    def key(a, b):
        i, j = sorted((index[a], index[b]))
        return i * len(names) + j

    # pairs whose values the panels draw from other tables: main c (the 200 isolates) and scale a (4,000 pairs)
    conventional = {g: pd.read_csv(f"{MTB}/results/conventional_tool/{g}.tsv", sep="\t") for g in GROUPS}
    conventional_keys = {key(a, b) for a, b in zip(conventional["all_cds"].A, conventional["all_cds"].B)}
    anchored = pd.read_csv(f"{trends}/containment_sample.tsv", sep="\t")
    anchored["key"] = [key(a, b) for a, b in zip(anchored.A, anchored.B)]
    shown = divergence_drawn[divergence_drawn.pairs >= parameters["divergence_min_pairs"]].reset_index(drop=True)
    divergence_bins = [tuple(shown.loc[k, ["divergence_low", "divergence_high"]]) for k in (0, len(shown) // 2, len(shown) - 1)]

    raw = {}
    for group in GROUPS:
        for scale in SCALES:
            print(f"reading {group}, scale factor {scale}", flush=True)
            wanted = (conventional_keys | set(anchored.key)) if scale == 1 else set()
            raw[group, scale] = raw_run(group, scale, args.jobs, names, wanted,
                                        divergence_bins if (group, scale) == ("all_cds", 1) else ())

    expected_pairs = len(names) * (len(names) - 1) // 2
    # ---- main a: the whole coding genome over every pair
    whole = raw["all_cds", 1]
    check("main a", "pairs drawn = n(n-1)/2 over the assembly-derived isolates", expected_pairs, whole["totals"]["pairs"], 0)
    check("main a", "caption: number of pairs", as_number(stated["mtbPairs"]), whole["totals"]["pairs"], 0)
    check("main a", "caption: percentage of pairs below 1", as_number(stated["mtbBelowOne"]),
          round(100 * whole["totals"]["below_one"] / whole["totals"]["finite"], 2), 0.005)
    histogram = pd.read_csv(f"{MTB}/results/tool/all_cds/scaled1/histogram.tsv", sep="\t")
    columns = ["same_sublineage", "same_lineage", "different_lineage"]
    check("main a", "finite values in the histogram drawn", int(histogram[columns].to_numpy().sum()), whole["totals"]["finite"], 0)
    inner = histogram.iloc[1:-1]
    centres = 0.5 * (inner.bin_low + inner.bin_high)
    outside = inner[(centres < 0.2) | (centres > 1.05)][columns].to_numpy().sum() / inner[columns].to_numpy().sum()
    checks.append(dict(panel="main a", check="share of pairs beyond the x axis (0.2 to 1.05), not drawn",
                       stated="", recomputed=outside, difference=0, tolerance=0, verdict="note"))

    # ---- main b: the boxes, from the exact quantiles
    for group in GROUPS:
        for name in QUANTILES:
            check("main b", f"{group} {name} drawn", groups.loc[group, name], raw[group, 1]["quantiles"][name], EXACT)
        check("main b", f"{group}: pairs", groups.loc[group, "pairs"], raw[group, 1]["totals"]["pairs"], 0)
    for group, macro in (("essential", "mtbMedianEssential"), ("all_cds", "mtbMedianWhole"),
                         ("nonessential", "mtbMedianNonessential"), ("pe_ppe", "mtbMedianPeppe"),
                         ("pe_ppe_hq", "mtbMedianPeppehq")):
        check("main b", f"text: median of {group} (3 decimals)", as_number(stated[macro]),
              round(raw[group, 1]["quantiles"]["median"], 3), 0.0005)
    for group, value in (("essential", 0.69), ("nonessential", 0.84), ("pe_ppe_hq", 0.90)):
        check("main b", f"diamond for {group} = Gan et al.'s per-gene pN/pS (Results, PMC12542625)",
              groups.loc[group, "published_per_gene"], value, EXACT)

    # ---- main c: against seqinr
    for group in GROUPS:
        table = conventional[group]
        table["key"] = [key(a, b) for a, b in zip(table.A, table.B)]
        rows = raw[group, 1]["wanted"].set_index("key")
        merged = table.join(rows, on="key", rsuffix="_raw")
        check("main c", f"{group}: pairs of the {len(set(table.A) | set(table.B))} isolates found in the raw output",
              len(table), merged.C_nt.notna().sum(), 0)
        # pairs identical over the group (both containments 1) have no value: 0/0 with either estimator
        undefined, undefined_raw = ~np.isfinite(merged.fmh_dnds), ~np.isfinite(merged.dnds)
        check("main c", f"{group}: pairs without a value in the table that have none in the raw output either",
              int(undefined.sum()), int((undefined & undefined_raw).sum()), 0)
        check("main c", f"{group}: pairs without a value in the raw output that have one in the table", 0,
              int((undefined_raw & ~undefined).sum()), 0)
        check("main c", f"{group}: pairs without a FracMinHash value that have a seqinr value", 0,
              int((undefined & np.isfinite(merged.seqinr_dnds)).sum()), 0)
        both = merged[["fmh_dnds", "dnds"]].replace([np.inf, -np.inf], np.nan).dropna()
        check("main c", f"{group}: largest |FracMinHash dN/dS in the table - raw|", 0, (both.fmh_dnds - both.dnds).abs().max(), EXACT)
        drawn = table.dropna(subset=["seqinr_dnds", "fmh_dnds"])
        check("main c", f"{group}: seqinr dN/dS = ka/ks", 0, (drawn.seqinr_dnds - drawn.ka / drawn.ks).abs().max(), EXACT)
    drawn = conventional["all_cds"].dropna(subset=["seqinr_dnds", "fmh_dnds"])
    limit = parameters["seqinr_limit"]
    beyond = ((drawn.seqinr_dnds > limit) | (drawn.fmh_dnds > limit) | (drawn.fmh_dnds < 0)).sum()
    check("main c", "caption: whole-genome pairs beyond the axes", as_number(stated["mtbSeqinrPanelOutside"]), beyond, 0)
    check("main c", "caption: pairs", as_number(stated["mtbSeqinrPairs"]), len(conventional["all_cds"]), 0)
    beyond_seqinr = (drawn.seqinr_dnds < 0).sum()
    check("main c", "whole-genome pairs with a negative seqinr value (would also be off the axes)", 0, beyond_seqinr, 0)
    random.seed(1)
    tasks = []
    for group in GROUPS:
        table = conventional[group].dropna(subset=["ka", "ks"])
        for _, row in table.sample(12, random_state=random.randrange(10**6)).iterrows():
            tasks.append((group, row.A, row.B, OUT))
    with Pool(min(len(tasks), args.jobs)) as pool:
        rerun = pd.DataFrame(pool.map(rerun_seqinr, tasks))
    rerun = rerun.merge(pd.concat([conventional[g].assign(group=g) for g in GROUPS])[["group", "A", "B", "ka", "ks"]],
                        on=["group", "A", "B"], suffixes=("_rerun", "_table"))
    rerun.to_csv(f"{OUT}/seqinr_rerun.tsv", sep="\t", index=False)
    check("main c", f"seqinr re-run on {len(rerun)} two-sequence alignments: largest |ka/ks difference|", 0,
          (rerun.ka_rerun / rerun.ks_rerun - rerun.ka_table / rerun.ks_table).abs().max(), 1e-6)

    # ---- main d: divergence bins
    for low, high in divergence_bins:
        values = whole["by_divergence"][(low, high)]
        row = shown[(shown.divergence_low == low) & (shown.divergence_high == high)].iloc[0]
        check("main d", f"bin {100 * low:.3f}-{100 * high:.3f}%: pairs", row.pairs, len(values), 0)
        for name, q in (("q1", 0.25), ("median", 0.5), ("q3", 0.75)):
            check("main d", f"bin {100 * low:.3f}-{100 * high:.3f}%: {name} drawn", row[name], np.quantile(values, q), EXACT)
    check("main d", "every bin drawn holds at least the stated number of pairs",
          as_number(stated["mtbDivergenceMinPairs"]), min(as_number(stated["mtbDivergenceMinPairs"]), shown.pairs.min()), 0)

    # ---- scale a: draft assemblies and anchored pairs
    rows = whole["wanted"].set_index("key")
    joined = anchored.join(rows, on="key", rsuffix="_raw")
    check("scale a", "anchored points: largest |C_nt drawn - raw|", 0, (joined.C_nt - joined.C_nt_raw).abs().max(), EXACT)
    check("scale a", "anchored points: largest |C_aa drawn - raw|", 0, (joined.C_aa - joined.C_aa_raw).abs().max(), EXACT)
    check("scale a", "caption: anchored points drawn", as_number(stated["mtbAnchoredPoints"]), joined.C_nt_raw.notna().sum(), 0)
    smoke = pd.read_csv(f"{MTB}/results/smoke/pairs.tsv", sep="\t")
    tool = pd.read_csv(f"{MTB}/work/smoke/scaled_1/fmh_omega_7.csv")
    tool["pair"] = [",".join(sorted((a, b))) for a, b in zip(tool.query_name, tool.match_name)]
    smoke["pair_sorted"] = [",".join(sorted((a, b))) for a, b in zip(smoke.A, smoke.B)]
    smoke_joined = smoke.merge(tool, left_on="pair_sorted", right_on="pair", suffixes=("", "_tool"))
    check("scale a", "draft-assembly pairs found in the tool's own output", len(smoke), len(smoke_joined), 0)
    check("scale a", "draft assemblies: largest |C_nt drawn - tool|", 0, (smoke_joined.C_nt_1 - smoke_joined.DNA_max_Cfrac).abs().max(), EXACT)
    check("scale a", "draft assemblies: largest |C_aa drawn - tool|", 0, (smoke_joined.C_aa_1 - smoke_joined.AA_max_Cfrac).abs().max(), EXACT)
    isolates = len(set(smoke.A) | set(smoke.B))
    check("scale a", "caption: draft-assembly pairs = n(n-1)/2", as_number(stated["mtbSmokePairs"]), isolates * (isolates - 1) // 2, 0)

    # ---- scale b, c, d: by scale factor, on this review's own subsample
    for group in GROUPS:
        exact = raw[group, 1]["sample"]
        for scale in SCALES:
            sketched = raw[group, scale]["sample"]
            both = pd.concat([exact, sketched], axis=1, keys=["exact", "sketched"], join="inner")
            both = both.replace([np.inf, -np.inf], np.nan).dropna()
            if scale != 1:
                check("scale b", f"{group}, scale {scale}: rank correlation with every k-mer",
                      scales.loc[(group, scale), "spearman_with_scale_1"], stats.spearmanr(both.exact, both.sketched)[0], SAMPLED)
            check("scale c", f"{group}, scale {scale}: median drawn", scales.loc[(group, scale), "median"],
                  raw[group, scale]["quantiles"]["median"], EXACT)
    for _, row in ordering.iterrows():
        lower = raw[row.lower, int(row.scale)]["sample"]
        higher = raw[row.higher, int(row.scale)]["sample"]
        both = pd.concat([lower, higher], axis=1, keys=["lower", "higher"], join="inner")
        both = both.replace([np.inf, -np.inf], np.nan).dropna()
        check("scale d", f"{row.lower} < {row.higher}, scale {int(row.scale)}: share of pairs", row.share_holding,
              (both.lower < both.higher).mean(), SAMPLED)

    result = pd.DataFrame(checks)
    result.to_csv(f"{OUT}/checks.tsv", sep="\t", index=False)
    counts = result.verdict.value_counts().to_dict()
    json.dump(counts, open(f"{OUT}/verdicts.json", "w"), indent=1)
    print(f"wrote {OUT}/checks.tsv: {counts}")


if __name__ == "__main__":
    main()
