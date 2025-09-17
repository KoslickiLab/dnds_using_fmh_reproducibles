#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import time
import subprocess
from pathlib import Path

def fmt_mmss(seconds_float: float) -> str:
    secs = int(round(seconds_float))
    return f"{secs // 60}:{secs % 60}"

def run_kaks_timed(kaks_bin: str, axt_in: Path, kaks_out: Path, method: str, log_path: Path):
    kaks_out.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [kaks_bin, "-i", str(axt_in), "-o", str(kaks_out), "-m", method]

    t0 = time.monotonic()
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         universal_newlines=True)
    dt = time.monotonic() - t0

    # write a simple log per run
    with open(log_path, "w") as lf:
        lf.write(f"$ {' '.join(cmd)}\n\n")
        lf.write("STDOUT:\n")
        lf.write(res.stdout or "")
        lf.write("\nSTDERR:\n")
        lf.write(res.stderr or "")

    if res.returncode != 0:
        raise RuntimeError(f"KaKs_Calculator failed (exit {res.returncode}). See log: {log_path}")

    return dt

def main():
    ap = argparse.ArgumentParser(description="Time KaKs_Calculator across sample sizes and methods.")
    ap.add_argument("--base_dir", default="/data/jzr5814/dnds_scalability_comparison")
    ap.add_argument("--samples", default="5,10,100,1000")
    ap.add_argument("--methods", default="NG,GY")
    ap.add_argument("--kaks_bin", default="KaKs_Calculator")
    ap.add_argument("--input_name", default="axt/sequences.axt")
    ap.add_argument("--out_template", default="axt/sequences_{method}.axt.kaks")
    ap.add_argument("--log_template", default="axt/sequences_{method}.kaks.log")
    ap.add_argument("--summary_name", default="timing_kaks_grid.tsv")
    args = ap.parse_args()

    base = Path(args.base_dir)
    base.mkdir(parents=True, exist_ok=True)

    samples = [s.strip() for s in args.samples.split(",") if s.strip()]
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]

    summary_path = base / args.summary_name
    with open(summary_path, "w") as sf:
        sf.write("\t".join(["sample", "method", "input_file", "output_file", "seconds", "status", "note"]) + "\n")

    total_ok = 0
    total_time = 0.0

    for s in samples:
        sample_dir = base / f"sample_{s}"
        in_path = sample_dir / args.input_name
        if not in_path.exists():
            note = f"missing input: {in_path}"
            print(f"[SKIP] {note}")
            with open(summary_path, "a") as sf:
                for m in methods:
                    sf.write("\t".join([s, m, str(in_path), "", "", "SKIP", note]) + "\n")
            continue

        for m in methods:
            out_path = sample_dir / args.out_template.format(method=m)
            log_path = sample_dir / args.log_template.format(method=m)
            try:
                secs = run_kaks_timed(args.kaks_bin, in_path, out_path, m, log_path)
                total_ok += 1
                total_time += secs
                print(f"[OK] sample_{s} {m}: {secs:.3f}s -> {out_path.name}")
                print(f"Mission accomplished. (Time elapsed: {fmt_mmss(secs)})")
                with open(summary_path, "a") as sf:
                    sf.write("\t".join([s, m, str(in_path), str(out_path), f"{secs:.6f}", "OK", ""]) + "\n")
            except Exception as e:
                note = str(e).replace("\n", " ")
                print(f"[FAIL] sample_{s} {m}: {note}")
                with open(summary_path, "a") as sf:
                    sf.write("\t".join([s, m, str(in_path), str(out_path), "", "FAIL", note]) + "\n")

    avg = (total_time / total_ok) if total_ok else 0.0
    print("\n==== KaKs Grid Timing ====")
    print(f"Completed: {total_ok} runs, total {total_time:.3f}s, avg {avg:.3f}s")
    print(f"Mission accomplished. (Time elapsed: {fmt_mmss(total_time)})")
    print(f"Summary: {summary_path.resolve()}")

if __name__ == "__main__":
    main()
