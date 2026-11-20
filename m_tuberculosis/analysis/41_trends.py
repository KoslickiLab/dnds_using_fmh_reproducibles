"""The comparisons of the vignette, from the summaries of the tool's output (step 40).

    usage:  python analysis/41_trends.py

  * groups.tsv        each group of genes at scale factor 1, every pair, beside Gan et al.'s values
  * paired.tsv        each pair of groups, pair by pair on the shared subsample, with a Wilcoxon test
  * scales.tsv        each group at every scale factor: median, quartiles, share above 1, the rank
                      correlation with scale factor 1 over the same pairs, and how many nucleotide
                      and protein hashes differ between the two isolates of a pair (median), which
                      is what the estimate rests on
  * ordering.tsv      at every scale factor, how often the orderings Gan et al. report hold pair by pair
  * whole_genome.tsv, divergence.tsv, sublineages.tsv   the whole coding genome at scale factor 1
  * containment_sample.tsv  4,000 pairs drawn at random from the subsample, with their containments,
                      for the containment panel of the figures
  * every_isolate.tsv the whole coding genome over the assembly- and read-derived isolates together
  * timings.tsv       the tool's own run times
  * prototype.tsv     the tool's distributions against those of the earlier prototype, bin by bin
  * summary.txt       all of the above as text
out: results/trends_tool/
"""
import csv
import glob
import json
import os
import zipfile

import numpy as np
import pandas as pd
from scipy import stats
from sourmash.minhash import _get_max_hash_for_scaled
from sourmash.signature import load_signatures_from_json

MTB = os.getcwd()
OUT = f"{MTB}/results/trends_tool"
os.makedirs(OUT, exist_ok=True)
GROUPS = ["essential", "all_cds", "nonessential", "pe_ppe", "pe_ppe_hq"]
SCALES = [1, 2, 3, 5, 10, 100, 1000]
# per-gene pN/pS reported by Gan et al. 2025; their PE/PPE value is over the high-quality segments
PUBLISHED = {"essential": 0.69, "nonessential": 0.84, "pe_ppe_hq": 0.90}
ORDERINGS = [("essential", "nonessential"), ("nonessential", "pe_ppe_hq"), ("essential", "pe_ppe_hq")]


def folder(run, scale):
    return f"{MTB}/results/tool/{run}/scaled{scale}"


def summary(run, scale):
    """Step 40's summary, with its quantiles (read from a histogram) replaced by the exact ones of
    step 40b."""
    path = f"{folder(run, scale)}/summary.tsv"
    if not os.path.exists(path):
        return None
    table = pd.read_csv(path, sep="\t").set_index("category")
    exact = pd.read_csv(f"{folder(run, scale)}/quantiles.tsv", sep="\t").set_index("category")
    assert (exact.pairs == table.pairs).all() and (exact.usable == table.usable).all()
    for column in ["p5", "q1", "median", "q3", "p95"]:
        table[column] = exact[column]
    return table


def subsample(run, scale):
    path = f"{folder(run, scale)}/subsample.parquet"
    return pd.read_parquet(path).set_index(["i", "j"]) if os.path.exists(path) else None


available = [g for g in GROUPS if summary(g, 1) is not None]
lines = []

# ---- groups at scale factor 1
rows = []
for group in available:
    row = summary(group, 1).loc["all pairs"]
    run = json.load(open(f"{folder(group, 1)}/run.json"))
    rows.append(dict(group=group, pairs=int(row.pairs), p5=row.p5, q1=row.q1, median=row["median"], q3=row.q3,
                     p95=row.p95, mean=row["mean"], below_one=row.below_one, above_one=row.above_one,
                     published_per_gene=PUBLISHED.get(group, np.nan), sketches=run["sketches"]))
groups = pd.DataFrame(rows)
groups.to_csv(f"{OUT}/groups.tsv", sep="\t", index=False)
lines += ["Groups of genes, scale factor 1 (every k-mer), every pair of the assembly-derived isolates",
          groups.round(4).to_string(index=False), ""]

# ---- pair by pair at scale factor 1
values = {group: subsample(group, 1).dnds for group in available}
rows = []
for a_index, a in enumerate(available):
    for b in available[a_index + 1:]:
        both = pd.concat([values[a], values[b]], axis=1, keys=[a, b], join="inner")
        both = both.replace([np.inf, -np.inf], np.nan).dropna()
        difference = both[b] - both[a]
        test_sample = difference.sample(min(len(difference), 200000), random_state=0)
        rows.append(dict(group_a=a, group_b=b, pairs=len(both), median_a=both[a].median(),
                         median_b=both[b].median(), median_difference=difference.median(),
                         share_b_greater=float((difference > 0).mean()),
                         wilcoxon_pairs=len(test_sample), wilcoxon_p=stats.wilcoxon(test_sample).pvalue))
paired = pd.DataFrame(rows)
paired.to_csv(f"{OUT}/paired.tsv", sep="\t", index=False)
lines += ["Pair by pair, scale factor 1, on the shared subsample", paired.round(5).to_string(index=False), ""]

# ---- scale factors
def reference_hashes(group):
    """The reference genome's sketch hashes (DNA, protein) of a group, read from the tool's sketches."""
    run = "all_cds_every_isolate" if group == "all_cds" else group
    out = {}
    for path in sorted(glob.glob(f"{MTB}/work/tool/{run}/sketch_shards/*/data.zip")):
        with zipfile.ZipFile(path) as archive:
            lines = [l for l in archive.read("SOURMASH-MANIFEST.csv").decode().splitlines() if not l.startswith("#")]
            for row in csv.DictReader(lines):
                if row["name"] == "ancestral_reference":
                    signature = next(iter(load_signatures_from_json(archive.read(row["internal_location"]))))
                    out[row["moltype"]] = np.fromiter(signature.minhash.hashes.keys(), dtype=np.uint64)
        if len(out) == 2:
            return out
    raise ValueError(f"no reference sketch found for {group}")


rows = []
for group in available:
    exact = subsample(group, 1)
    reference = reference_hashes(group)
    for scale in SCALES:
        table = summary(group, scale)
        if table is None:
            continue
        row = table.loc["all pairs"]
        sketched = subsample(group, scale)
        both = pd.concat([exact.dnds, sketched.dnds], axis=1, keys=["exact", "sketched"], join="inner")
        both = both.replace([np.inf, -np.inf], np.nan).dropna()
        # how many hashes of a pair differ, after sketching: what the estimate is computed from
        max_hash = np.uint64(_get_max_hash_for_scaled(scale))
        protein_size = int((reference["protein"] <= max_hash).sum())
        dna_size = int((reference["DNA"] <= max_hash).sum())
        rows.append(dict(group=group, scale=scale, median=row["median"], q1=row.q1, q3=row.q3,
                         p5=row.p5, p95=row.p95, above_one=row.above_one,
                         usable=row.usable / row.pairs,
                         spearman_with_scale_1=both.exact.corr(both.sketched, method="spearman")
                         if scale != 1 else 1.0,
                         median_absolute_error=(both.sketched - both.exact).abs().median(),
                         reference_protein_hashes=protein_size, reference_dna_hashes=dna_size,
                         differing_protein_hashes=float(((1 - sketched.C_aa) * protein_size).median()),
                         differing_dna_hashes=float(((1 - sketched.C_nt) * dna_size).median())))
scales = pd.DataFrame(rows)
scales.to_csv(f"{OUT}/scales.tsv", sep="\t", index=False)
lines += ["Every group at every scale factor", scales.round(4).to_string(index=False), ""]

rows = []
for scale in SCALES:
    for a, b in ORDERINGS:
        sa, sb = subsample(a, scale), subsample(b, scale)
        if sa is None or sb is None:
            continue
        both = pd.concat([sa.dnds, sb.dnds], axis=1, keys=["a", "b"], join="inner")
        both = both.replace([np.inf, -np.inf], np.nan).dropna()
        rows.append(dict(scale=scale, lower=a, higher=b, pairs=len(both),
                         share_holding=float((both.b > both.a).mean()),
                         median_lower=both.a.median(), median_higher=both.b.median()))
ordering = pd.DataFrame(rows)
ordering.to_csv(f"{OUT}/ordering.tsv", sep="\t", index=False)
lines += ["How often each published ordering holds pair by pair, at each scale factor",
          ordering.round(4).to_string(index=False), ""]

# ---- whole coding genome
whole = summary("all_cds", 1)
whole.reset_index().to_csv(f"{OUT}/whole_genome.tsv", sep="\t", index=False)
names = [line.strip() for line in open(f"{folder('all_cds', 1)}/names.txt") if line.strip()]
sample = pd.read_parquet(f"{folder('all_cds', 1)}/subsample.parquet").sample(4000, random_state=0)
sample = sample.assign(A=[names[i] for i in sample.i], B=[names[j] for j in sample.j])
sample[["A", "B", "C_nt", "C_aa", "dnds"]].to_csv(f"{OUT}/containment_sample.tsv", sep="\t", index=False)
divergence = pd.read_csv(f"{folder('all_cds', 1)}/quantiles_by_divergence.tsv", sep="\t")
divergence.to_csv(f"{OUT}/divergence.tsv", sep="\t", index=False)
by_sublineage = pd.read_csv(f"{folder('all_cds', 1)}/by_sublineage.tsv", sep="\t")
within = by_sublineage[by_sublineage.sublineage_a == by_sublineage.sublineage_b].rename(
    columns={"sublineage_a": "sublineage"})[["sublineage", "pairs", "mean"]].sort_values("mean")
within.to_csv(f"{OUT}/sublineages.tsv", sep="\t", index=False)
lines += ["Whole coding genome by category of pair, scale factor 1", whole.round(4).to_string(), "",
          "Whole coding genome against per-nucleotide divergence", divergence.round(5).to_string(index=False), "",
          "Within-sublineage means", within.round(4).to_string(index=False), ""]

every = summary("all_cds_every_isolate", 1)
if every is not None:
    every.reset_index().to_csv(f"{OUT}/every_isolate.tsv", sep="\t", index=False)
    lines += ["Whole coding genome, assembly- and read-derived isolates together", every.round(4).to_string(), ""]

# ---- run times, from the tool's run.json
rows = []
for run in ["all_cds_every_isolate"] + available:
    for scale in SCALES:
        path = f"{folder(run, scale)}/run.json"
        if os.path.exists(path):
            report = json.load(open(path))
            rows.append(dict(run=run, scale=scale, sketches=report["sketches"], pairs=report["pairs"],
                             cores=report["cores"], load_seconds=report["load_seconds"],
                             pairs_seconds=report["pairs_seconds"],
                             core_seconds_per_million_pairs=report["cores"] * report["pairs_seconds"]
                             / report["pairs"] * 1e6))
timings = pd.DataFrame(rows)
timings.to_csv(f"{OUT}/timings.tsv", sep="\t", index=False)
lines += ["The tool's run times", timings.round(3).to_string(index=False), ""]

# ---- against the earlier prototype (same pairs, same bins: the counts should be the same)
rows = []
for group in available:
    for scale in SCALES:
        prototype_path = f"{MTB}/results/prototype/all_pairs/{group}" + ("" if scale == 1 else f"_scale{scale}") + "/histogram.tsv"
        tool_path = f"{folder(group, scale)}/histogram.tsv"
        if not (os.path.exists(prototype_path) and os.path.exists(tool_path)):
            continue
        prototype = pd.read_csv(prototype_path, sep="\t")
        tool = pd.read_csv(tool_path, sep="\t")
        columns = ["same_sublineage", "same_lineage", "different_lineage"]
        difference = (prototype[columns] - tool[columns]).abs()
        rows.append(dict(group=group, scale=scale, prototype_pairs=int(prototype[columns].to_numpy().sum()),
                         tool_pairs=int(tool[columns].to_numpy().sum()),
                         pairs_in_different_bins=int(difference.to_numpy().sum() // 2)))
prototype = pd.DataFrame(rows)
prototype.to_csv(f"{OUT}/prototype.tsv", sep="\t", index=False)
lines += ["The tool against the earlier prototype, bin by bin (0.005 wide)", prototype.to_string(index=False)]

open(f"{OUT}/summary.txt", "w").write("\n".join(lines) + "\n")
pd.set_option("display.width", 250)
print("\n".join(lines))
