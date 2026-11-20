"""The runs whose reads go through the read pipeline (steps 07 and 08).

    usage:  python analysis/06_read_runs.py [out_folder]

Two kinds of run: every run of Gan et al. without an AllTheBacteria assembly, and a calibration set
of runs that do have one (12 per sublineage, drawn at random among the runs that ENA lists), on which
the substitutions found from reads and from the assembly can be compared (step 34).

out (data/isolates/ by default):
    calibration_runs.tsv     run, sublineage and ATB sample of the 300 calibration runs
    read_pipeline_runs.txt   the calibration runs, then the runs without an assembly
"""
import os
import sys

import pandas as pd

OUT = sys.argv[1] if len(sys.argv) > 1 else "data/isolates"
os.makedirs(OUT, exist_ok=True)
matched = pd.read_csv("data/isolates/gan_to_atb.tsv", sep="\t", dtype=str)
ena = pd.read_csv("data/isolates/ena_runinfo.tsv", sep="\t", dtype=str)

with_assembly = matched[matched["file"].notna() & matched.run.isin(ena.run_accession)]
calibration = with_assembly.groupby("sublineage").sample(n=12, random_state=7)
calibration[["run", "sublineage", "sample"]].to_csv(f"{OUT}/calibration_runs.tsv", sep="\t", index=False)
without_assembly = matched[matched["file"].isna()].run.drop_duplicates()
runs = list(calibration.run) + list(without_assembly)
open(f"{OUT}/read_pipeline_runs.txt", "w").write("\n".join(runs) + "\n")
print(len(calibration), "calibration runs;", len(without_assembly), "runs without an assembly;", len(runs), "in all")
