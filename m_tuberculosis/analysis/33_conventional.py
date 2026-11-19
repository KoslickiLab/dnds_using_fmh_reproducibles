"""The conventional estimate Gan et al. 2025 used (R package seqinr, kaks) against FracMinHash dN/dS,
on exactly the same sequences.

    usage:  python analysis/33_conventional.py [n_isolates] [seed]

For n assembly-derived isolates chosen at random (200 and seed 20261118 by default, i.e. 19,900
pairs; seqinr needs about 0.05 s per pair of whole coding genomes on one core, so this takes a
quarter of an hour), and for each group of genes:
  * seqinr's kaks is run on the coding sequences of step 30 concatenated per isolate. Every isolate
    has the same coordinates, so the concatenations are already aligned. seqinr's own run time is
    measured inside R;
  * the FracMinHash dN/dS of the same pairs is read from the tool's run over all pairs (step 32,
    scale factor 1), so both estimates come from the same sequences.

out: results/conventional_tool/isolates.txt        the isolates
     results/conventional_tool/<group>.tsv         per pair: seqinr ka, ks, ka/ks and the tool's values
     results/conventional_tool/summary.tsv         per group: medians, correlations, seqinr's time
"""
import glob
import os
import random
import subprocess
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

MTB = os.getcwd()
GROUPS = ["all_cds", "essential", "nonessential", "pe_ppe", "pe_ppe_hq"]
N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 20261118
WORK = f"{MTB}/work/conventional_tool"
OUT = f"{MTB}/results/conventional_tool"
os.makedirs(WORK, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

assembly = sorted(line.strip() for line in open(f"{MTB}/work/fasta/assembly_names.txt") if line.strip())
random.seed(SEED)
names = sorted(random.sample(assembly, N))
open(f"{OUT}/isolates.txt", "w").write("\n".join(names) + "\n")

R_SCRIPT = """
suppressMessages(library(seqinr))
args <- commandArgs(trailingOnly=TRUE)
alignment <- read.alignment(args[1], format="fasta")
start <- Sys.time()
result <- kaks(alignment)
seconds <- as.numeric(difftime(Sys.time(), start, units="secs"))
ka <- as.matrix(result$ka); ks <- as.matrix(result$ks)
pairs <- which(upper.tri(ka), arr.ind=TRUE)
out <- data.frame(A=alignment$nam[pairs[,1]], B=alignment$nam[pairs[,2]], ka=ka[pairs], ks=ks[pairs])
write.table(out, args[2], sep="\\t", row.names=FALSE, quote=FALSE)
cat(sprintf("%.3f\\n", seconds), file=args[3])
"""
open(f"{WORK}/kaks.R", "w").write(R_SCRIPT)


def concatenated(path):
    return "".join(line.strip() for line in open(path) if not line.startswith(">"))


def pairs_in(task):
    """The tool's rows for pairs of two wanted isolates, from one parquet part."""
    path, wanted = task
    table = pq.read_table(path, columns=["query_name", "match_name", "DNA_max_Cfrac", "AA_max_Cfrac",
                                         "PdN", "PdS", "dN/dS"]).to_pandas()
    query, match = table.query_name.astype(str), table.match_name.astype(str)
    return table[query.isin(wanted) & match.isin(wanted)]


if __name__ == "__main__":
    # seqinr, all groups side by side (one core each)
    running = {}
    for group in GROUPS:
        with open(f"{WORK}/{group}.fasta", "w") as out:
            for name in names:
                out.write(f">{name}\n{concatenated(f'{MTB}/work/fasta/{group}/{name}.ffn')}\n")
        running[group] = subprocess.Popen(["Rscript", f"{WORK}/kaks.R", f"{WORK}/{group}.fasta",
                                           f"{WORK}/{group}.kaks.tsv", f"{WORK}/{group}.seconds"],
                                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(f"seqinr running on {N} isolates ({N * (N - 1) // 2} pairs) for {len(GROUPS)} groups", flush=True)

    # the tool's values for the same pairs
    wanted = set(names)
    tool = {}
    for group in GROUPS:
        parts = sorted(glob.glob(f"{MTB}/work/tool/{group}/reference_scaled1/pairs/part-*.parquet"))
        with Pool(100) as pool:
            found = pd.concat(pool.map(pairs_in, [(p, wanted) for p in parts]), ignore_index=True)
        found["pair"] = [" | ".join(sorted((a, b)))
                         for a, b in zip(found.query_name.astype(str), found.match_name.astype(str))]
        tool[group] = found.set_index("pair")
        print(f"  {group}: {len(found)} pairs from the tool's run", flush=True)

    rows = []
    for group in GROUPS:
        process = running[group]
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            sys.exit(f"seqinr failed for {group}: {stderr[-1500:]}")
        conventional = pd.read_csv(f"{WORK}/{group}.kaks.tsv", sep="\t")
        conventional["pair"] = [" | ".join(sorted((a, b))) for a, b in zip(conventional.A, conventional.B)]
        conventional["seqinr_dnds"] = conventional.ka / conventional.ks
        both = conventional.set_index("pair").join(
            tool[group][["DNA_max_Cfrac", "AA_max_Cfrac", "PdN", "PdS", "dN/dS"]], how="inner")
        both = both.rename(columns={"dN/dS": "fmh_dnds"}).replace([np.inf, -np.inf], np.nan)
        both.reset_index().to_csv(f"{OUT}/{group}.tsv", sep="\t", index=False)
        usable = both.dropna(subset=["seqinr_dnds", "fmh_dnds"])
        seconds = float(open(f"{WORK}/{group}.seconds").read())
        rows.append(dict(group=group, isolates=N, seed=SEED, pairs=len(both), usable=len(usable),
                         median_seqinr=usable.seqinr_dnds.median(), median_fmh=usable.fmh_dnds.median(),
                         spearman=usable.seqinr_dnds.corr(usable.fmh_dnds, method="spearman"),
                         pearson=usable.seqinr_dnds.corr(usable.fmh_dnds),
                         seqinr_seconds=seconds, seqinr_seconds_per_pair=seconds / len(conventional)))
    summary = pd.DataFrame(rows)
    summary.to_csv(f"{OUT}/summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(summary.round(4).to_string(index=False))
