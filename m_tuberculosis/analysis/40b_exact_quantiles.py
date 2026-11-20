"""Exact quantiles of the FracMinHash dN/dS values of every run and scale factor.

    usage:  python analysis/40b_exact_quantiles.py [--jobs N] [run ...]

Step 40 reduces each run to a histogram with bins 0.005 wide, which places a quantile only to within
half a bin. This step reads the tool's per-pair values again and computes the quantiles exactly
(numpy's default: linear interpolation between order statistics), over the finite values, with the
categories and divergence bins of step 40. Pairs that involve the reference genome are left out.

out (results/tool/<run>/scaled<S>/):
    quantiles.tsv                 per category and for all pairs: pairs, usable, p5, q1, median, q3, p95
    quantiles_by_divergence.tsv   whole coding genome only: per bin of divergence, 1 - C_nt^(1/21),
                                  the number of usable pairs and their q1, median and q3
"""
import argparse
import glob
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

MTB = os.getcwd()
REFERENCE_NAME = "ancestral_reference"
CATEGORIES = ["same sublineage", "same lineage, different sublineage", "different lineage"]
DIVERGENCE_BINS = np.concatenate([np.linspace(0, 0.001, 51), [np.inf]])     # as in step 40
QUANTILES = {"p5": 0.05, "q1": 0.25, "median": 0.5, "q3": 0.75, "p95": 0.95}
RUNS = ["all_cds", "essential", "nonessential", "pe_ppe", "pe_ppe_hq", "all_cds_every_isolate"]
WITH_DIVERGENCE = {"all_cds"}


def initialise(names, sublineage, lineage):
    global NAME_INDEX, SUBLINEAGE, LINEAGE
    pa.set_cpu_count(1)              # one thread per worker process
    pa.set_io_thread_count(1)
    NAME_INDEX = {name: index for index, name in enumerate(names)}
    SUBLINEAGE, LINEAGE = sublineage, lineage


def positions(column):
    """Positions in names.txt of a dictionary-encoded name column."""
    parts = []
    for chunk in column.chunks:
        mapping = np.array([NAME_INDEX[n] for n in chunk.dictionary.to_pylist()], dtype=np.int64)
        parts.append(mapping[chunk.indices.to_numpy()])
    return np.concatenate(parts)


def one_part(path):
    """The finite dN/dS values of one part, with the category and divergence bin of each pair, and
    the number of pairs per category (finite or not)."""
    table = pq.read_table(path, columns=["query_name", "match_name", "DNA_max_Cfrac", "dN/dS"], use_threads=False)
    query, match = positions(table.column("query_name")), positions(table.column("match_name"))
    dnds = table.column("dN/dS").to_numpy()
    c_nt = table.column("DNA_max_Cfrac").to_numpy()
    isolates = (SUBLINEAGE[query] >= 0) & (SUBLINEAGE[match] >= 0)
    category = np.where(SUBLINEAGE[query] == SUBLINEAGE[match], 0,
                        np.where(LINEAGE[query] == LINEAGE[match], 1, 2)).astype(np.int8)
    pairs = np.bincount(category[isolates], minlength=3)
    usable = isolates & np.isfinite(dnds)
    divergence = 1.0 - np.power(np.clip(c_nt[usable], 0, None), 1.0 / 21)
    divergence_bin = np.clip(np.digitize(divergence, DIVERGENCE_BINS) - 1, 0, len(DIVERGENCE_BINS) - 2)
    return pairs, dnds[usable], category[usable], divergence_bin.astype(np.int8)


def labels(names):
    """Sublineage and lineage index of every name (-1 for the reference), from Gan et al.'s Table S1
    as matched to the runs and samples of AllTheBacteria (step 05)."""
    matched = pd.read_csv(f"{MTB}/data/isolates/gan_to_atb.tsv", sep="\t", dtype=str)
    by_name = dict(zip(matched.run, matched.sublineage))
    by_name.update({s: l for s, l in zip(matched["sample"], matched.sublineage) if isinstance(s, str)})
    sublineages = sorted({by_name[n] for n in names if n != REFERENCE_NAME})
    lineages = sorted({s[:2] for s in sublineages})
    sublineage = np.array([sublineages.index(by_name[n]) if n != REFERENCE_NAME else -1 for n in names])
    lineage = np.array([lineages.index(by_name[n][:2]) if n != REFERENCE_NAME else -1 for n in names])
    return sublineage, lineage


def exact(values):
    if len(values) == 0:
        return {name: np.nan for name in QUANTILES}
    return dict(zip(QUANTILES, np.quantile(values, list(QUANTILES.values()))))


def run_one(run, folder, jobs):
    scale = folder.rsplit("reference_scaled", 1)[1]
    out = f"{MTB}/results/tool/{run}/scaled{scale}"
    names = [line.strip() for line in open(f"{folder}/names.txt") if line.strip()]
    sublineage, lineage = labels(names)
    parts = sorted(glob.glob(f"{folder}/pairs/part-*.parquet"))
    pairs = np.zeros(3, dtype=np.int64)
    values, categories, divergence_bins = [], [], []
    with Pool(jobs, initializer=initialise, initargs=(names, sublineage, lineage)) as pool:
        for p, v, c, d in pool.imap_unordered(one_part, parts, chunksize=1):
            pairs += p
            values.append(v)
            categories.append(c)
            divergence_bins.append(d)
    values, categories = np.concatenate(values), np.concatenate(categories)

    rows = [dict(category=name, pairs=int(pairs[index]), usable=int((categories == index).sum()),
                 **exact(values[categories == index])) for index, name in enumerate(CATEGORIES)]
    rows.append(dict(category="all pairs", pairs=int(pairs.sum()), usable=len(values), **exact(values)))
    pd.DataFrame(rows).to_csv(f"{out}/quantiles.tsv", sep="\t", index=False)

    if run in WITH_DIVERGENCE:
        divergence_bins = np.concatenate(divergence_bins)
        rows = []
        for index in range(len(DIVERGENCE_BINS) - 1):
            inside = values[divergence_bins == index]
            row = exact(inside)
            rows.append(dict(divergence_low=DIVERGENCE_BINS[index], divergence_high=DIVERGENCE_BINS[index + 1],
                             pairs=len(inside), q1=row["q1"], median=row["median"], q3=row["q3"]))
        pd.DataFrame(rows).to_csv(f"{out}/quantiles_by_divergence.tsv", sep="\t", index=False)
    print(f"{run} scale {scale}: {len(values):,} usable values, median {np.median(values):.4f}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("runs", nargs="*", default=RUNS)
    parser.add_argument("--jobs", type=int, default=200, help="worker processes, one thread each")
    parser.add_argument("--scales", type=int, nargs="*", help="only these scale factors (default: all)")
    args = parser.parse_args()
    for run in args.runs:
        for folder in sorted(glob.glob(f"{MTB}/work/tool/{run}/reference_scaled*")):
            scale = int(folder.rsplit("reference_scaled", 1)[1])
            if os.path.isdir(f"{folder}/pairs") and (not args.scales or scale in args.scales):
                run_one(run, folder, args.jobs)


if __name__ == "__main__":
    main()
