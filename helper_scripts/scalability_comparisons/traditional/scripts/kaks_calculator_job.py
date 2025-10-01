#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import time
import subprocess
from pathlib import Path

def fmt_mmss(seconds_float: float) -> str:
    secs = int(round(seconds_float))
    return f"{secs // 60}:{secs % 60}"

def run_kaks_timed(kaks_bin: str, axt_in: Path, kaks_out: Path, method: str, log_path: Path) -> float:
    kaks_out.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [kaks_bin, "-i", str(axt_in), "-o", str(kaks_out), "-m", method]

    t0 = time.monotonic()
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    dt = time.monotonic() - t0

    with open(log_path, "w") as lf:
        lf.write(f"$ {' '.join(cmd)}\n\nSTDOUT:\n{res.stdout or ''}\nSTDERR:\n{res.stderr or ''}")

    if res.returncode != 0:
        raise RuntimeError(f"KaKs_Calculator failed (exit {res.returncode}). See log: {log_path}")

    return dt

def main():
    ap = argparse.ArgumentParser(description="Time KaKs_Calculator across sample sizes and methods; optionally for every individual AXT.")
    ap.add_argument("--base_dir", default="/data/jzr5814/dnds_scalability_comparison", help="Base directory containing sample_<N>/")
    ap.add_argument("--samples", default="5,10,100,1000", help="Comma-separated sample sizes (numbers only)")
    ap.add_argument("--methods", default="NG,GY", help="Comma-separated KaKs methods (e.g. NG,GY,MS,MA)")
    ap.add_argument("--kaks_bin", default="KaKs_Calculator", help="KaKs_Calculator binary path or name")

    # Single combined AXT (existing behavior)
    ap.add_argument("--input_name", default="axt/sequences.axt", help="Input AXT path relative to each sample dir")
    ap.add_argument("--out_template", default="axt/sequences_{method}.axt.kaks", help="(filename or path) template per sample")
    ap.add_argument("--log_template", default="axt/sequences_{method}.kaks.log", help="(filename or path) template per sample")
    ap.add_argument("--summary_name", default="timing_kaks_grid.tsv", help="Summary TSV at base_dir for single combined runs")
    ap.add_argument("--skip_single", action="store_true", help="Skip the single combined AXT run")

    # Per-AXT batch timing
    ap.add_argument("--axt_glob", default=None, help="Glob under each sample dir to find MANY AXT files (e.g. 'axt/*.axt' or 'axt/**/*.axt')")
    ap.add_argument("--summary_name_files", default="timing_kaks_per_axt.tsv", help="Summary TSV at base_dir for per-AXT runs")

    # NEW: put outputs/logs into a 'kaks' subdir
    ap.add_argument("--kaks_subdir", default="kaks", help="Subfolder (relative to sample or AXT dir) to store .kaks + logs")

    args = ap.parse_args()

    base = Path(args.base_dir)
    base.mkdir(parents=True, exist_ok=True)

    samples = [s.strip() for s in args.samples.split(",") if s.strip()]
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]

    # Summaries
    grid_summary = base / args.summary_name
    with open(grid_summary, "w") as sf:
        sf.write("\t".join(["sample", "method", "input_file", "output_file", "seconds", "status", "note"]) + "\n")

    per_file_summary = base / args.summary_name_files
    with open(per_file_summary, "w") as sf:
        sf.write("\t".join(["sample", "method", "axt_file", "kaks_file", "seconds", "status", "note"]) + "\n")

    total_ok = 0
    total_time = 0.0
    total_file_ok = 0
    total_file_time = 0.0

    for s in samples:
        sample_dir = base / f"sample_{s}"

        # --------- Single combined AXT run (optional) ----------
        if not args.skip_single:
            in_path = sample_dir / args.input_name
            if not in_path.exists():
                note = f"missing input: {in_path}"
                print(f"[SKIP] {note}")
                with open(grid_summary, "a") as sf:
                    for m in methods:
                        sf.write("\t".join([s, m, str(in_path), "", "", "SKIP", note]) + "\n")
            else:
                for m in methods:
                    # Write into sample_<N>/<kaks_subdir>/, using only the basename of the templates
                    out_name = Path(args.out_template.format(method=m)).name
                    log_name = Path(args.log_template.format(method=m)).name
                    out_path = sample_dir / args.kaks_subdir / out_name
                    log_path = sample_dir / args.kaks_subdir / log_name

                    try:
                        secs = run_kaks_timed(args.kaks_bin, in_path, out_path, m, log_path)
                        total_ok += 1
                        total_time += secs
                        print(f"[OK] sample_{s} {m} (combined): {secs:.3f}s -> {out_path}")
                        print(f"Mission accomplished. (Time elapsed: {fmt_mmss(secs)})")
                        with open(grid_summary, "a") as sf:
                            sf.write("\t".join([s, m, str(in_path), str(out_path), f"{secs:.6f}", "OK", ""]) + "\n")
                    except Exception as e:
                        note = str(e).replace("\n", " ")
                        print(f"[FAIL] sample_{s} {m} (combined): {note}")
                        with open(grid_summary, "a") as sf:
                            sf.write("\t".join([s, m, str(in_path), str(out_path), "", "FAIL", note]) + "\n")

        # --------- Per-AXT runs (each individual AXT) ----------
        if args.axt_glob:
            matches = list(sample_dir.glob(args.axt_glob))
            if not matches:
                print(f"[SKIP] No AXT matched in sample_{s} with pattern '{args.axt_glob}'")
            else:
                for axt_in in sorted(matches):
                    if not axt_in.is_file() or axt_in.suffix.lower() != ".axt":
                        continue
                    # outputs/logs go into <axt_in.parent>/<kaks_subdir>/
                    out_dir = axt_in.parent / args.kaks_subdir
                    for m in methods:
                        kaks_name = f"{axt_in.stem}_{m}.axt.kaks"
                        log_name  = f"{axt_in.stem}_{m}.kaks.log"
                        kaks_out = out_dir / kaks_name
                        log_path = out_dir / log_name
                        try:
                            secs = run_kaks_timed(args.kaks_bin, axt_in, kaks_out, m, log_path)
                            total_file_ok += 1
                            total_file_time += secs
                            print(f"[OK] sample_{s} {m} (per-file): {secs:.3f}s -> {kaks_out}")
                            with open(per_file_summary, "a") as sf:
                                sf.write("\t".join([s, m, str(axt_in), str(kaks_out), f"{secs:.6f}", "OK", ""]) + "\n")
                        except Exception as e:
                            note = str(e).replace("\n", " ")
                            print(f"[FAIL] sample_{s} {m} (per-file) {axt_in.name}: {note}")
                            with open(per_file_summary, "a") as sf:
                                sf.write("\t".join([s, m, str(axt_in), str(kaks_out), "", "FAIL", note]) + "\n")

    # Totals
    avg = (total_time / total_ok) if total_ok else 0.0
    print("\n==== KaKs Grid Timing (combined files) ====")
    print(f"Completed: {total_ok} runs, total {total_time:.3f}s, avg {avg:.3f}s")
    print(f"Mission accomplished. (Time elapsed: {fmt_mmss(total_time)})")
    print(f"Summary: {grid_summary.resolve()}")

    avg_file = (total_file_time / total_file_ok) if total_file_ok else 0.0
    print("\n==== KaKs Timing (per AXT files) ====")
    print(f"Completed: {total_file_ok} runs, total {total_file_time:.3f}s, avg {avg_file:.3f}s")
    print(f"Summary: {per_file_summary.resolve()}")

if __name__ == "__main__":
    main()
