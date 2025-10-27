#!/usr/bin/env python3
import os
import subprocess
import argparse
import pandas as pd
import time

def run_command(cmd):
    start_time = time.time()
    subprocess.run(cmd, shell=True, check=True)
    return time.time() - start_time

def main():
    parser = argparse.ArgumentParser(description="Compare DNA and protein signatures with sourmash.")
    parser.add_argument("--base", required=True, help="Base directory with sample_X/rep_Y folders.")
    parser.add_argument("--ksize", type=int, default=7, help="K-mer size (default: 7)")
    parser.add_argument("--outcsv", required=True, help="Output CSV for all timing records")
    args = parser.parse_args()

    samples = [5, 10]
    replicates = range(1, 11)

    all_records = []

    for s in samples:
        for r in replicates:
            base_dir = os.path.join(args.base, f"sample_{s}", f"rep_{r}")
            sig_dir = os.path.join(base_dir, "signatures")
            dna_ref = os.path.join(sig_dir, f"sample{s}_rep{r}_dna.sig")
            prot_ref = os.path.join(sig_dir, f"sample{s}_rep{r}_protein.sig")

            for molecule, ref in [("dna", dna_ref), ("protein", prot_ref)]:
                out_csv = os.path.join(base_dir, "containments", f"compare.{molecule}.{args.ksize}.csv")
                os.makedirs(os.path.dirname(out_csv), exist_ok=True)

                cmd = f"sourmash compare {ref} {ref} --containment --{molecule} --ksize {args.ksize} --csv {out_csv}"
                duration = run_command(cmd)

                all_records.append({
                    "sample": s,
                    "replicate": r,
                    "molecule": molecule,
                    "ksize": args.ksize,
                    "duration_sec": round(duration, 2)
                })

    df = pd.DataFrame(all_records)
    df.to_csv(args.outcsv, index=False)
    print(f"✅ All comparison timings saved to {args.outcsv}")

if __name__ == "__main__":
    main()
