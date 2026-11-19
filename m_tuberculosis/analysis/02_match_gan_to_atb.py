"""Match the 49,954 runs of Gan et al. 2025 (Table S1) to AllTheBacteria 0.2 assemblies.

    usage:  python analysis/02_match_gan_to_atb.py
    out:  data/isolates/gan_runs.tsv           Table S1 as a plain table (run, sublineage, BioProject)
          data/isolates/gan_to_atb.tsv         one row per run: ATB sample, ATB species call, file, HQ, stats, n runs of the sample
          results/data_overlap/summary.txt
          results/data_overlap/collection.tsv  distinct runs of Table S1, and those with an ATB assembly
"""
import pandas as pd, gzip, os
MTB = os.getcwd()
os.makedirs("data/isolates", exist_ok=True)
s1 = pd.read_excel("data/gan2025/supp/files/msystems.00898-25-s0008.xlsx", header=1)
s1.columns = ["run", "sublineage", "bioproject"]
s1.to_csv("data/isolates/gan_runs.tsv", sep="\t", index=False)
meta = "data/atb/0.2/metadata/"
ena = pd.read_csv(meta + "ena_metadata.tsv.gz", sep="\t", usecols=["sample_accession", "run_accession"], dtype=str)
runs_per_sample = ena.groupby("sample_accession").size().rename("runs_of_sample")
s2s = pd.read_csv(meta + "sample2species2file.tsv.gz", sep="\t", dtype=str)
calls = pd.read_csv(meta + "species_calls.tsv.gz", sep="\t", dtype=str).rename(columns={"Sample": "sample", "Species": "species_call", "HQ": "hq"})
stats = pd.read_csv(meta + "assembly-stats.tsv.gz", sep="\t", usecols=["sample", "total_length", "number", "N50"], dtype={"sample": str})
ck = pd.read_csv(meta + "checkm2.tsv.gz", sep="\t", usecols=["sample", "Completeness_Specific", "Contamination"], dtype={"sample": str})
m = s1.merge(ena.rename(columns={"run_accession": "run", "sample_accession": "sample"}), on="run", how="left")
m = m.merge(runs_per_sample, left_on="sample", right_index=True, how="left")
m = m.merge(s2s, on="sample", how="left").merge(calls, on="sample", how="left").merge(stats, on="sample", how="left").merge(ck, on="sample", how="left")
m.to_csv("data/isolates/gan_to_atb.tsv", sep="\t", index=False)
lines = [f"Gan et al. Table S1 runs: {len(s1)} (unique {s1.run.nunique()})",
         f"runs found in ATB ENA metadata: {m['sample'].notna().sum()}",
         f"runs with an ATB 0.2 assembly: {m['file'].notna().sum()}",
         f"  of which ATB species call M. tuberculosis: {(m.species_call == 'Mycobacterium tuberculosis').sum()}",
         f"  of which in the ATB high-quality set: {(m.hq == 'T').sum()}",
         f"  samples with >1 run in ENA (ATB assembles per sample): {(m.runs_of_sample > 1).sum()}",
         f"  distinct ATB samples: {m.loc[m.file.notna(), 'sample'].nunique()}",
         "species calls among matched runs:", m.species_call.value_counts(dropna=False).head(10).to_string(),
         "sublineage counts (all runs):", s1.sublineage.value_counts().to_string(),
         "sublineage counts (runs with ATB assembly):", m.loc[m.file.notna()].sublineage.value_counts().to_string()]
open("results/data_overlap/summary.txt", "w").write("\n".join(lines) + "\n")
# the two counts the vignette quotes (step 43)
pd.DataFrame([dict(runs=s1.run.nunique(), runs_with_assembly=int(m["file"].notna().sum()))]).to_csv(
    "results/data_overlap/collection.tsv", sep="\t", index=False)
print("\n".join(lines))
