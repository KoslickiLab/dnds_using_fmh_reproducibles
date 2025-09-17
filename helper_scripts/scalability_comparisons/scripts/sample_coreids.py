#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re, random, argparse
from pathlib import Path

CORE_ID_RE = re.compile(r'(cds_WP_\d+\.\d)')

def read_fasta(path):
    header, seq = None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq)
                header, seq = line.strip(), []
            else:
                seq.append(line.strip())
        if header is not None:
            yield header, ''.join(seq)

def wrap_seq(s, width=60):
    return '\n'.join(s[i:i+width] for i in range(0, len(s), width))

def core_id_from_header(h):
    m = CORE_ID_RE.search(h); return m.group(1) if m else None

def index_by_core_id(fa):
    idx = {}
    for h, s in read_fasta(fa):
        cid = core_id_from_header(h)
        if cid:
            idx[cid] = (h, s)
    return idx

def main():
    ap = argparse.ArgumentParser(description="Sample common cds_WP IDs across 3 strains and write 3-seq FASTAs.")
    ap.add_argument("--strain1", required=True)
    ap.add_argument("--strain2", required=True)
    ap.add_argument("--strain3", required=True)
    ap.add_argument("--sample_n", type=int, required=True, help="How many gene sets to sample")
    ap.add_argument("--min_len", type=int, default=1000, help="Length filter on strain1 only")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out_dir", required=True, help="Destination folder (e.g., .../sample_N/rep_R/fasta)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    idx1 = index_by_core_id(args.strain1)
    idx2 = index_by_core_id(args.strain2)
    idx3 = index_by_core_id(args.strain3)

    eligible_in_1 = {cid for cid,(h,s) in idx1.items() if len(s) > args.min_len}
    common = eligible_in_1 & idx2.keys() & idx3.keys()
    if not common:
        raise SystemExit(f"No common CDS (>{args.min_len}) across three strains.")

    ids = sorted(common)
    if args.sample_n > len(ids):
        print(f"[WARN] requested {args.sample_n}, only {len(ids)} available; sampling all.")
        args.sample_n = len(ids)

    random.seed(args.seed)
    picked = random.sample(ids, args.sample_n)

    manifest = ["core_id\toutfile\tstrain_file\theader"]
    for cid in picked:
        out_fa = out_dir / f"{cid}.fasta"
        with open(out_fa, "w") as out:
            for src, idx in [(args.strain1, idx1),(args.strain2, idx2),(args.strain3, idx3)]:
                h, s = idx[cid]
                if not h.startswith(">"): h = f">{h}"
                out.write(f"{h}\n{wrap_seq(s)}\n")
                manifest.append(f"{cid}\t{out_fa.name}\t{Path(src).name}\t{h[1:]}")

    (out_dir / "manifest.tsv").write_text("\n".join(manifest) + "\n")
    (out_dir / "_SAMPLED.OK").write_text("\n")    # marker for Snakemake
    print(f"[OK] wrote {len(picked)} FASTAs to {out_dir}")

if __name__ == "__main__":
    main()
