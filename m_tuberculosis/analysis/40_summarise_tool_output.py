"""Summaries of the FracMinHash dN/dS tool's output (step 32), one per run and scale factor.

    usage:  python analysis/40_summarise_tool_output.py [--jobs N] [run ...]

Reads work/tool/<run>/reference_scaled<S>/pairs/*.parquet, the per-pair table the tool writes (the
columns of its main mode's fmh_omega_7.csv), and reduces it to the distributions the vignette uses.
Nothing is recomputed: dN/dS, PdN, PdS and the containments are the tool's values. Pairs that involve
the reference genome itself are left out; it is not one of the isolates.

Pairs are put in three categories by the sublineages Gan et al. assigned (Table S1): both isolates in
one sublineage; one lineage, different sublineages; different lineages. Divergence is the per-
nucleotide estimate of the manuscript, 1 - C_nt^(1/21).

out (results/tool/<run>/scaled<S>/):
    summary.tsv        per category: pairs, usable values, mean, quantiles, share below 1, share above 1
    histogram.tsv      the distribution of dN/dS per category (same bins as the earlier prototype)
    by_divergence.tsv  quartiles of dN/dS per bin of divergence
    by_sublineage.tsv  mean dN/dS per pair of sublineages
    subsample.parquet  every pair (i, j) with (i * 2654435761 + j * 97) mod 64 == 0: the same pairs in
                       every run, so that runs can be compared pair by pair
    run.json           the tool's own run report, copied
"""
import argparse
import glob
import json
import os
import shutil
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

MTB = os.getcwd()
REFERENCE_NAME = "ancestral_reference"
BINS = np.concatenate([[-np.inf], np.linspace(-2, 8, 2001), [np.inf]])
DIVERGENCE_BINS = np.concatenate([np.linspace(0, 0.001, 51), [np.inf]])
SUBSAMPLE_MODULUS = 64
CATEGORIES = ["same sublineage", "same lineage, different sublineage", "different lineage"]


def labels_for(names):
    """Sublineage index and lineage index of every name (-1 for the reference)."""
    matched = pd.read_csv(f"{MTB}/data/isolates/gan_to_atb.tsv", sep="\t", dtype=str)
    by_name = {**dict(zip(matched.run, matched.sublineage)),
               **{s: l for s, l in zip(matched["sample"], matched.sublineage) if isinstance(s, str)}}
    sublineages = sorted({by_name[n] for n in names if n != REFERENCE_NAME})
    lineages = sorted({s[:2] for s in sublineages})
    sublineage = np.array([sublineages.index(by_name[n]) if n != REFERENCE_NAME else -1 for n in names])
    lineage = np.array([lineages.index(by_name[n][:2]) if n != REFERENCE_NAME else -1 for n in names])
    return sublineages, sublineage, lineage


def global_indices(column):
    """Positions in the tool's names.txt of a dictionary-encoded name column (each file, and each chunk
    of it, may carry its own dictionary)."""
    parts = []
    for chunk in column.chunks:
        mapping = np.array([NAME_INDEX[n] for n in chunk.dictionary.to_pylist()], dtype=np.int64)
        parts.append(mapping[chunk.indices.to_numpy()])
    return np.concatenate(parts)


def one_part(path):
    table = pq.read_table(path, columns=["query_name", "match_name", "DNA_max_Cfrac", "AA_max_Cfrac", "dN/dS"])
    query = global_indices(table.column("query_name"))
    match = global_indices(table.column("match_name"))
    c_nt = table.column("DNA_max_Cfrac").to_numpy()
    c_aa = table.column("AA_max_Cfrac").to_numpy()
    dnds = table.column("dN/dS").to_numpy()
    keep = (SUBLINEAGE[query] >= 0) & (SUBLINEAGE[match] >= 0)
    query, match, c_nt, c_aa, dnds = query[keep], match[keep], c_nt[keep], c_aa[keep], dnds[keep]

    same_sublineage = SUBLINEAGE[query] == SUBLINEAGE[match]
    same_lineage = LINEAGE[query] == LINEAGE[match]
    category = np.where(same_sublineage, 0, np.where(same_lineage, 1, 2))
    finite = np.isfinite(dnds)
    histogram = np.zeros((3, len(BINS) - 1), dtype=np.int64)
    totals = np.zeros((3, 6))
    for which in range(3):
        selected = category == which
        usable = selected & finite
        histogram[which] = np.histogram(dnds[selected & finite], bins=BINS)[0]
        values = dnds[usable]
        totals[which] = [selected.sum(), usable.sum(), values.sum(), (values ** 2).sum(),
                         (values < 1).sum(), (values > 1).sum()]

    divergence = 1.0 - np.power(np.clip(c_nt, 0, None), 1.0 / 21)
    d_index = np.clip(np.digitize(divergence, DIVERGENCE_BINS) - 1, 0, len(DIVERGENCE_BINS) - 2)
    v_index = np.clip(np.digitize(dnds, BINS) - 1, 0, len(BINS) - 2)
    by_divergence = np.bincount(d_index[finite] * (len(BINS) - 1) + v_index[finite],
                                minlength=(len(DIVERGENCE_BINS) - 1) * (len(BINS) - 1))

    n_sub = SUBLINEAGE.max() + 1
    low = np.minimum(SUBLINEAGE[query], SUBLINEAGE[match])
    high = np.maximum(SUBLINEAGE[query], SUBLINEAGE[match])
    cell = (low * n_sub + high)[finite]
    by_sublineage = np.stack([np.bincount(cell, minlength=n_sub * n_sub),
                              np.bincount(cell, weights=dnds[finite], minlength=n_sub * n_sub)])

    picked = (query * 2654435761 + match * 97) % SUBSAMPLE_MODULUS == 0
    sample = pd.DataFrame(dict(i=query[picked].astype(np.int32), j=match[picked].astype(np.int32),
                               C_nt=c_nt[picked], C_aa=c_aa[picked], dnds=dnds[picked],
                               category=category[picked].astype(np.int8)))
    return histogram, totals, by_divergence, by_sublineage, sample


def initialise(sublineage, lineage, names):
    global SUBLINEAGE, LINEAGE, NAME_INDEX
    SUBLINEAGE, LINEAGE = sublineage, lineage
    NAME_INDEX = {name: index for index, name in enumerate(names)}


def quantiles(counts, centres, qs):
    cumulative = np.cumsum(counts)
    if cumulative[-1] == 0:
        return [np.nan] * len(qs)
    return [centres[np.searchsorted(cumulative, q * cumulative[-1])] for q in qs]


def summarise(run, folder, jobs):
    scale = folder.rsplit("reference_scaled", 1)[1]
    out = f"{MTB}/results/tool/{run}/scaled{scale}"
    os.makedirs(out, exist_ok=True)
    parts = sorted(glob.glob(f"{folder}/pairs/part-*.parquet"))
    names = [line.strip() for line in open(f"{folder}/names.txt") if line.strip()]
    sublineages, sublineage, lineage = labels_for(names)
    histogram = np.zeros((3, len(BINS) - 1), dtype=np.int64)
    totals = np.zeros((3, 6))
    by_divergence = np.zeros((len(DIVERGENCE_BINS) - 1) * (len(BINS) - 1), dtype=np.int64)
    by_sublineage = np.zeros((2, len(sublineages) ** 2))
    samples = []
    with Pool(jobs, initializer=initialise, initargs=(sublineage, lineage, names)) as pool:
        for h, t, d, s, sample in pool.imap_unordered(one_part, parts, chunksize=1):
            histogram += h
            totals += t
            by_divergence += d
            by_sublineage += s
            samples.append(sample)

    centres = 0.5 * (BINS[1:-2] + BINS[2:-1])
    rows = []
    for which, name in enumerate(CATEGORIES):
        count, usable, total, total_squared, below, above = totals[which]
        p5, q1, median, q3, p95 = quantiles(histogram[which][1:-1], centres, (0.05, 0.25, 0.5, 0.75, 0.95))
        rows.append(dict(category=name, pairs=int(count), usable=int(usable),
                         mean=total / usable if usable else np.nan,
                         p5=p5, q1=q1, median=median, q3=q3, p95=p95,
                         below_one=below / usable if usable else np.nan,
                         above_one=above / usable if usable else np.nan))
    all_counts = histogram.sum(axis=0)[1:-1]
    p5, q1, median, q3, p95 = quantiles(all_counts, centres, (0.05, 0.25, 0.5, 0.75, 0.95))
    count, usable = totals[:, 0].sum(), totals[:, 1].sum()
    rows.append(dict(category="all pairs", pairs=int(count), usable=int(usable),
                     mean=totals[:, 2].sum() / usable, p5=p5, q1=q1, median=median, q3=q3, p95=p95,
                     below_one=totals[:, 4].sum() / usable, above_one=totals[:, 5].sum() / usable))
    pd.DataFrame(rows).to_csv(f"{out}/summary.tsv", sep="\t", index=False)
    pd.DataFrame(dict(bin_low=BINS[:-1], bin_high=BINS[1:], same_sublineage=histogram[0],
                      same_lineage=histogram[1], different_lineage=histogram[2])).to_csv(
        f"{out}/histogram.tsv", sep="\t", index=False)
    grid = by_divergence.reshape(len(DIVERGENCE_BINS) - 1, len(BINS) - 1)
    divergence_rows = []
    for d in range(len(DIVERGENCE_BINS) - 1):
        counts = grid[d][1:-1]
        if counts.sum():
            q1, median, q3 = quantiles(counts, centres, (0.25, 0.5, 0.75))
            divergence_rows.append(dict(divergence_low=DIVERGENCE_BINS[d], divergence_high=DIVERGENCE_BINS[d + 1],
                                        pairs=int(counts.sum()), q1=q1, median=median, q3=q3))
    pd.DataFrame(divergence_rows).to_csv(f"{out}/by_divergence.tsv", sep="\t", index=False)
    n = len(sublineages)
    sublineage_rows = [dict(sublineage_a=sublineages[a], sublineage_b=sublineages[b],
                            pairs=int(by_sublineage[0, a * n + b]),
                            mean=by_sublineage[1, a * n + b] / by_sublineage[0, a * n + b])
                       for a in range(n) for b in range(a, n) if by_sublineage[0, a * n + b]]
    pd.DataFrame(sublineage_rows).to_csv(f"{out}/by_sublineage.tsv", sep="\t", index=False)
    pd.concat(samples, ignore_index=True).to_parquet(f"{out}/subsample.parquet", index=False)
    shutil.copy(f"{folder}/names.txt", f"{out}/names.txt")
    shutil.copy(f"{folder}/run.json", f"{out}/run.json")
    print(f"{run} scale {scale}: {int(count)} pairs, median {rows[-1]['median']:.4f}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs", type=int, default=100)
    parser.add_argument("runs", nargs="*")
    arguments = parser.parse_args()
    runs = arguments.runs or sorted(os.listdir(f"{MTB}/work/tool"))
    for run in runs:
        for folder in sorted(glob.glob(f"{MTB}/work/tool/{run}/reference_scaled*"),
                             key=lambda f: int(f.rsplit("reference_scaled", 1)[1])):
            if glob.glob(f"{folder}/pairs/part-*.parquet"):
                summarise(run, folder, arguments.jobs)
