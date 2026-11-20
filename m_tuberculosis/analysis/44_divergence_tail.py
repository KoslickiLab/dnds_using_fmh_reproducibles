"""Who makes up the most divergent pairs of the whole-coding-genome comparison?

    usage:  python analysis/44_divergence_tail.py [threshold]

Pairs of isolates whose per-nucleotide divergence (1 - C_nt^(1/21)) exceeds the threshold (default
0.0006, i.e. 0.06%) are rare and carry the highest values in the divergence panel. This lists the
isolates that take part in them, on the subsample of step 40, with each isolate's number of
substitutions against the reference, and repeats the divergence summary without those isolates.

out: results/trends_tool/divergence_tail.tsv          isolates in the tail, most frequent first
     results/trends_tool/divergence_tail_summary.tsv  the numbers the vignette quotes
     results/trends_tool/divergence_without_tail.tsv  quartiles by divergence without them
"""
import os
import sys

import numpy as np
import pandas as pd

MTB = os.getcwd()
THRESHOLD = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0006
folder = f"{MTB}/results/tool/all_cds/scaled1"
names = np.array([line.strip() for line in open(f"{folder}/names.txt") if line.strip()])
sample = pd.read_parquet(f"{folder}/subsample.parquet")
sample["divergence"] = 1 - np.power(sample.C_nt.clip(lower=0), 1 / 21)
inputs = pd.read_csv(f"{MTB}/results/inputs/isolates.tsv", sep="\t").set_index("name")
counts = pd.Series(np.bincount(np.concatenate([sample.i, sample.j]), minlength=len(names)))

tail = sample[sample.divergence > THRESHOLD]
involved = pd.Series(np.concatenate([tail.i, tail.j])).value_counts()
table = pd.DataFrame({"isolate": names[involved.index], "pairs_in_tail": involved.values,
                      "pairs_in_subsample": counts[involved.index].values})
table["share_of_its_pairs_in_tail"] = table.pairs_in_tail / table.pairs_in_subsample
table["substitutions"] = table.isolate.map(inputs.substitutions)
outliers = table[table.share_of_its_pairs_in_tail > 0.5]
table.to_csv(f"{MTB}/results/trends_tool/divergence_tail.tsv", sep="\t", index=False)
print(f"{len(tail)} of {len(sample)} subsampled pairs ({100 * len(tail) / len(sample):.3f}%) above "
      f"{100 * THRESHOLD:.3f}% divergence; {len(table)} isolates take part, {len(outliers)} of them in the "
      f"tail for most of their pairs")
print(f"substitutions of those {len(outliers)} isolates: median {outliers.substitutions.median():.0f} "
      f"(all isolates: median {inputs.substitutions.median():.0f})")
outlier_index = set(involved.index[(table.share_of_its_pairs_in_tail > 0.5).to_numpy()])
in_tail_with_outlier = (tail.i.isin(outlier_index) | tail.j.isin(outlier_index)).mean()
print(f"share of tail pairs involving at least one of them: {100 * in_tail_with_outlier:.1f}%")
print(table.head(12).to_string(index=False))
pd.DataFrame([dict(threshold=THRESHOLD, tail_pairs=len(tail), subsample_pairs=len(sample),
                   tail_share=len(tail) / len(sample), outlier_isolates=len(outliers),
                   outlier_median_substitutions=outliers.substitutions.median(),
                   all_median_substitutions=inputs.substitutions.median(),
                   tail_share_with_outlier=in_tail_with_outlier)]).to_csv(
    f"{MTB}/results/trends_tool/divergence_tail_summary.tsv", sep="\t", index=False)

kept = sample[~(sample.i.isin(outlier_index) | sample.j.isin(outlier_index))]
bins = np.linspace(0, 0.001, 51)
kept = kept.assign(bin=np.digitize(kept.divergence, bins) - 1)
rows = [dict(divergence_low=bins[b], divergence_high=bins[b + 1], pairs=len(g),
             q1=g.dnds.quantile(0.25), median=g.dnds.median(), q3=g.dnds.quantile(0.75))
        for b, g in kept.groupby("bin") if 0 <= b < 50 and len(g) >= 10000 // 64]
pd.DataFrame(rows).to_csv(f"{MTB}/results/trends_tool/divergence_without_tail.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).round(5).to_string(index=False))
