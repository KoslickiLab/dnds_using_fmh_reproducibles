#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, csv, os
from pathlib import Path
from typing import Optional, Dict, List

def read_tsv_sum_seconds(path: Path, seconds_col: str) -> float:
    """Sum a float column from a headered TSV; returns 0.0 if file missing."""
    if not path.exists():
        return 0.0
    total = 0.0
    with path.open() as f:
        rdr = csv.DictReader(f, delimiter="\t")
        if seconds_col not in rdr.fieldnames:
            # try to accept 'seconds ' with whitespace etc.
            # map normalized keys
            cols = {c.strip().lower(): c for c in rdr.fieldnames or []}
            key = cols.get(seconds_col.strip().lower())
            if key is None:
                raise ValueError(f"{path} has no column '{seconds_col}' (found: {rdr.fieldnames})")
            for row in rdr:
                total += float(row[key] or 0.0)
            return total
        for row in rdr:
            total += float(row[seconds_col] or 0.0)
    return total

def sum_kaks_from_benchmarks(kaks_root: Path) -> float:
    """Sum 's' (seconds) from Snakemake benchmark TSVs under kaks_root."""
    if not kaks_root.exists():
        return 0.0
    total = 0.0
    found = False
    for bench in kaks_root.rglob("*.benchmark.txt"):
        with bench.open() as f:
            header = f.readline().strip().split("\t")
            row = f.readline().strip().split("\t")
            if not header or not row or len(header) != len(row):
                continue
            try:
                idx = header.index("s")
                total += float(row[idx])
                found = True
            except ValueError:
                # some environments write other headers; ignore
                continue
    return total if found else 0.0

def sum_kaks_from_global_tsv(global_tsv: Path, sample: str, rep: Optional[str]) -> float:
    """Fallback: sum seconds from timing_kaks_per_axt.tsv filtered by sample[/rep]."""
    if not global_tsv.exists():
        return 0.0
    total = 0.0
    with global_tsv.open() as f:
        rdr = csv.DictReader(f, delimiter="\t")
        # Expect columns: sample, method, axt_file, kaks_file, seconds, status, note
        for row in rdr:
            if row.get("status") != "OK":
                continue
            if row.get("sample") != str(sample):
                continue
            axt_file = row.get("axt_file", "")
            if rep and f"/rep_{rep}/" not in axt_file:
                continue
            try:
                total += float(row.get("seconds", 0.0))
            except ValueError:
                pass
    return total

def main():
    ap = argparse.ArgumentParser(description="Summarize timing across alignment, AXT conversion, and KaKs.")
    ap.add_argument("--base_dir", required=True, help="Root directory (contains sample_<N>/ …)")
    ap.add_argument("--samples", default="5,10,100,1000", help="Comma-separated sample sizes, e.g. 5,10,100,1000")
    ap.add_argument("--replicates", type=int, default=1, help="Replicates per sample (if using rep_<r> subfolders)")
    ap.add_argument("--out", default="timing_summary_overall.tsv", help="Output TSV path (relative or absolute)")
    ap.add_argument("--prefer_benchmarks", action="store_true",
                    help="Use Snakemake benchmark files for KaKs if present (recommended)")
    ap.add_argument("--global_kaks_tsv", default=None,
                    help="Optional path to timing_kaks_per_axt.tsv (fallback if no benchmarks)")
    args = ap.parse_args()

    base = Path(args.base_dir)
    samples = [s.strip() for s in args.samples.split(",") if s.strip()]
    out_path = Path(args.out)
    if not out_path.is_absolute():
        out_path = base / out_path

    # header
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as out:
        out.write("\t".join([
            "sample", "replicate",
            "align_seconds", "axt_seconds", "kaks_seconds",
            "total_seconds"
        ]) + "\n")

    for s in samples:
        for r in range(1, args.replicates + 1):
            # Build paths (with or without replicates)
            if args.replicates > 1:
                align_tsv = base / f"sample_{s}" / f"rep_{r}" / "alignments" / "timing_alignments.tsv"
                axt_tsv   = base / f"sample_{s}" / f"rep_{r}" / "axt" / "timing_summary.tsv"
                kaks_root = base / f"sample_{s}" / f"rep_{r}" / "kaks"
            else:
                align_tsv = base / f"sample_{s}" / "alignments" / "timing_alignments.tsv"
                axt_tsv   = base / f"sample_{s}" / "axt" / "timing_summary.tsv"
                kaks_root = base / f"sample_{s}" / "kaks"

            # 1) alignment time
            align_seconds = read_tsv_sum_seconds(align_tsv, "seconds")

            # 2) axt conversion time
            axt_seconds = read_tsv_sum_seconds(axt_tsv, "convert_seconds")

            # 3) KaKs time
            kaks_seconds = 0.0
            used_bench = False
            if args.prefer_benchmarks:
                kaks_seconds = sum_kaks_from_benchmarks(kaks_root)
                used_bench = kaks_seconds > 0.0
            if not used_bench:
                # fallback: global per-axt timing file
                if args.global_kaks_tsv:
                    rep = str(r) if args.replicates > 1 else None
                    kaks_seconds = sum_kaks_from_global_tsv(Path(args.global_kaks_tsv), str(s), rep)

            total = align_seconds + axt_seconds + kaks_seconds

            with out_path.open("a") as out:
                out.write("\t".join([
                    str(s), str(r),
                    f"{align_seconds:.6f}",
                    f"{axt_seconds:.6f}",
                    f"{kaks_seconds:.6f}",
                    f"{total:.6f}",
                ]) + "\n")

    print(f"[OK] wrote {out_path}")

if __name__ == "__main__":
    main()
