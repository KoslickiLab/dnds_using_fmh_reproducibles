"""Smoke test: FracMinHash dN/dS between M. tuberculosis isolates, at five scale factors.

    usage:  python analysis/12_smoke_scales.py [n_per_sublineage]

Runs the published implementation (external/dnds-using-fmh/fmh_dnds.py, unmodified) on a small set of
isolates, once per scale factor, with the parameters of the manuscript (k_aa = 7, so k_nt = 21).
Isolates are the highest-quality assembly of each of the 25 sublineages of Gan et al. Table S1.

The question this answers: isolates of one species differ at only ~0.1% of sites, so the estimator's
denominator (C_aa^(1/7) - C_nt^(1/7)) is a difference of two small numbers. Sketching adds noise to
both. This measures how much, by comparing each scale factor with scale factor 1, which keeps every
k-mer and is therefore exact.

out: results/smoke/scale_<s>.tsv      per-pair values at scale factor s
     results/smoke/summary.tsv        one row per scale factor
     results/smoke/pairs.tsv          per-pair values of every scale factor, side by side
"""
import os
import subprocess
import sys
import numpy as np
import pandas as pd

MTB = os.getcwd()
LAB = f"{MTB}/external/dnds-using-fmh"
SCALES = [1, 10, 100, 500, 1000]
N_PER_SUBLINEAGE = int(sys.argv[1]) if len(sys.argv) > 1 else 2
OUT = f"{MTB}/results/smoke"
WORK = f"{MTB}/work/smoke"
os.makedirs(OUT, exist_ok=True)
os.makedirs(WORK, exist_ok=True)

# ---- isolates: the best assembly of each sublineage, among those with gene calls
isolates = pd.read_csv(f"{MTB}/data/isolates/gan_to_atb.tsv", sep="\t", dtype=str)
isolates = isolates[(isolates.species_call == "Mycobacterium tuberculosis") & (isolates.hq == "T")]
isolates["n_contigs"] = pd.to_numeric(isolates.number)
isolates["contamination"] = pd.to_numeric(isolates.Contamination)
have = {f[:-7] for f in os.listdir(f"{MTB}/work/genes") if f.endswith(".faa.gz")}
isolates = isolates[isolates["sample"].isin(have)].drop_duplicates("sample")
chosen = (isolates.sort_values(["contamination", "n_contigs"])
                  .groupby("sublineage").head(N_PER_SUBLINEAGE)
                  .sort_values(["sublineage", "sample"]))
print(f"{len(chosen)} isolates from {chosen.sublineage.nunique()} sublineages")
chosen[["sample", "sublineage", "run", "n_contigs", "contamination"]].to_csv(f"{OUT}/isolates.tsv", sep="\t", index=False)

sheet = ["name,genome_filename,protein_filename"] + [
    f"{s},{MTB}/work/genes/{s}.ffn.gz,{MTB}/work/genes/{s}.faa.gz" for s in chosen["sample"]]
open(f"{WORK}/samples.csv", "w").write("\n".join(sheet) + "\n")

# ---- the published implementation, once per scale factor
sublineage = dict(zip(chosen["sample"], chosen.sublineage))
tables = {}
for scale in SCALES:
    directory = f"{WORK}/scaled_{scale}"
    os.makedirs(directory, exist_ok=True)
    for stale in ("data.zip", f"fmh_omega_7.csv"):
        if os.path.exists(f"{directory}/{stale}"):
            os.remove(f"{directory}/{stale}")
    command = [sys.executable, f"{LAB}/fmh_dnds.py", "--fasta_input_list", f"{WORK}/samples.csv",
               "--ksize", "7", "--scaled_input", str(scale), "--directory", directory,
               "--cores", os.environ.get("THREADS", "64"), "--threshold", "0"]
    log = open(f"{directory}/run.log", "w")
    elapsed = -1.0
    import time
    start = time.time()
    env = dict(os.environ)
    result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, cwd=LAB, env=env)
    elapsed = time.time() - start
    if result.returncode != 0:
        sys.exit(f"fmh_dnds.py failed at scale {scale}; see {directory}/run.log")
    table = pd.read_csv(f"{directory}/fmh_omega_7.csv")
    table = table.rename(columns={"query_name": "A", "match_name": "B", "DNA_max_Cfrac": "C_nt",
                                  "AA_max_Cfrac": "C_aa", "dN/dS": "fmh_dnds"})
    table = table[table.A != table.B].copy()
    table["pair"] = [",".join(sorted((a, b))) for a, b in zip(table.A, table.B)]
    table = table.drop_duplicates("pair")
    table["same_sublineage"] = [sublineage[a] == sublineage[b] for a, b in zip(table.A, table.B)]
    table["scale"] = scale
    table["seconds"] = elapsed
    table.to_csv(f"{OUT}/scale_{scale}.tsv", sep="\t", index=False)
    tables[scale] = table
    print(f"scale {scale:5d}: {len(table)} pairs in {elapsed:.0f} s")

# ---- how each scale factor compares with scale factor 1 (every k-mer, exact)
exact = tables[1].set_index("pair")
rows = []
side_by_side = exact[["A", "B", "same_sublineage"]].copy()
for scale, table in tables.items():
    joined = table.set_index("pair").reindex(exact.index)
    side_by_side[f"dnds_{scale}"] = joined.fmh_dnds
    side_by_side[f"C_nt_{scale}"] = joined.C_nt
    side_by_side[f"C_aa_{scale}"] = joined.C_aa
    finite = joined.fmh_dnds.replace([np.inf, -np.inf], np.nan)
    usable = finite.notna() & (finite > 0)
    error = (finite - exact.fmh_dnds) / exact.fmh_dnds
    rows.append(dict(scale=scale, pairs=len(table), seconds=table.seconds.iloc[0],
                     median_dnds=finite.median(), q1=finite.quantile(0.25), q3=finite.quantile(0.75),
                     fraction_usable=usable.mean(),
                     median_abs_relative_error=error.abs().median() if scale != 1 else 0.0,
                     spearman_with_exact=finite.corr(exact.fmh_dnds, method="spearman"),
                     median_within_sublineage=finite[exact.same_sublineage].median(),
                     median_between_sublineage=finite[~exact.same_sublineage].median(),
                     median_C_nt=joined.C_nt.median(), median_C_aa=joined.C_aa.median()))
summary = pd.DataFrame(rows)
summary.to_csv(f"{OUT}/summary.tsv", sep="\t", index=False)
side_by_side.to_csv(f"{OUT}/pairs.tsv", sep="\t")
pd.set_option("display.width", 200)
print(summary.to_string(index=False))
