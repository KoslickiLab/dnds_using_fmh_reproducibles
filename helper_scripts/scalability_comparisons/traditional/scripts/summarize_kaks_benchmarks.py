#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
import glob

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Glob pattern for all benchmark files")
parser.add_argument("--output", required=True, help="TSV summary output file")
args = parser.parse_args()

out_path = Path(args.output)
out_path.parent.mkdir(parents=True, exist_ok=True)

all_files = glob.glob(args.input_glob, recursive=True)

with open(out_path, "w", newline="") as csvfile:
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(["sample", "replicate", "gene", "pair", "method", "user_time_sec", "system_time_sec", "elapsed_time_sec"])

    for bench_file in sorted(all_files):
        parts = Path(bench_file).parts
        sample = [p for p in parts if p.startswith("sample_")][0].split("_")[1]
        rep = [p for p in parts if p.startswith("rep_")][0].split("_")[1]
        gene = Path(bench_file).parent.name
        stem = Path(bench_file).stem
        pair_method = stem.split("_")
        pair = pair_method[0]
        method = pair_method[1]

        user, sys, elapsed = None, None, None
        with open(bench_file) as f:
            for line in f:
                if "User time" in line:
                    user = float(line.split()[3])
                elif "System time" in line:
                    sys = float(line.split()[3])
                elif "Elapsed (wall clock) time" in line:
                    elapsed = float(line.split()[6])
        writer.writerow([sample, rep, gene, pair, method, user, sys, elapsed])
