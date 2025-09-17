#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import shutil
import subprocess
import time
from pathlib import Path

# Try clustalw2 first, then clustalw
def which_clustalw():
    for exe in ("clustalw2", "clustalw"):
        p = shutil.which(exe)
        if p:
            return p
    return None

def read_fasta(path):
    """Yield (header, sequence) for each record in FASTA."""
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

def sanitize_label(h):
    """
    Make a compact label from a header: prefer the first whitespace-delimited token
    without the initial '>' and strip odd chars.
    """
    token = h[1:].split()[0] if h.startswith('>') else h.split()[0]
    # keep alnum, dot, bar, underscore, dash
    return re.sub(r'[^A-Za-z0-9._|\-]+', '_', token)

def write_pair_fasta(out_path: Path, recA, recB):
    """Write a 2-seq FASTA file. rec = (header, seq)."""
    with open(out_path, "w") as out:
        for h, s in (recA, recB):
            h2 = h if h.startswith('>') else f">{h}"
            out.write(h2 + "\n")
            for i in range(0, len(s), 60):
                out.write(s[i:i+60] + "\n")

def run_clustalw(clustalw_bin: str, infile: Path, outfile: Path, is_dna=True, quiet=True):
    args = [
        clustalw_bin,
        f"-INFILE={infile}",
        "-ALIGN",
        f"-OUTFILE={outfile}",
        "-OUTPUT=CLUSTAL",
        "-OUTORDER=INPUT",
    ]
    if is_dna:
        args.append("-TYPE=DNA")
    if quiet:
        args.append("-QUIET")

    res = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         universal_newlines=True)
    if res.returncode != 0:
        raise RuntimeError(f"ClustalW failed on {infile}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")

    # Clean the .dnd neighbor-joining guide tree ClustalW leaves
    dnd = infile.with_suffix(".dnd")
    if dnd.exists():
        try: dnd.unlink()
        except Exception: pass

def enumerate_pairs(records, include_self=True, unique=False):
    """
    records: list[(header, seq)]
    - include_self=True  -> i,j over all indices (N*N)
    - include_self=False -> i != j (ordered) -> N*(N-1)
    - unique=True        -> unordered unique pairs (i<j); overrides include_self
    """
    n = len(records)
    if unique:
        for i in range(n):
            for j in range(i+1, n):
                yield i, j
        return
    for i in range(n):
        for j in range(n):
            if not include_self and i == j:
                continue
            yield i, j

def ensure_summary(path: Path):
    """Create summary TSV with header if it doesn't exist."""
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as sf:
            sf.write("\t".join([
                "fasta_file",
                "i_idx",
                "j_idx",
                "labelA",
                "labelB",
                "pair_fasta",
                "aln_file",
                "seconds",
                "status",
                "note"
            ]) + "\n")

def append_summary(path: Path, row: list):
    with open(path, "a") as sf:
        sf.write("\t".join(row) + "\n")

def fmt_mmss(seconds_float: float) -> str:
    secs = int(round(seconds_float))
    return f"{secs // 60}:{secs % 60}"

def main():
    ap = argparse.ArgumentParser(description="Run pairwise ClustalW alignments within each FASTA.")
    ap.add_argument("--input", required=True,
                    help="Path to a FASTA file OR a directory containing FASTA files.")
    ap.add_argument("--out_dir", required=True,
                    help="Directory to write alignments into.")
    ap.add_argument("--protein", action="store_true",
                    help="Treat sequences as protein (default DNA).")
    ap.add_argument("--no-self", action="store_true",
                    help="Exclude self-comparisons; gives N*(N-1) alignments.")
    ap.add_argument("--unique", action="store_true",
                    help="Use unordered unique pairs (i<j); gives N choose 2 alignments.")
    ap.add_argument("--exts", nargs="+", default=[".fa", ".fna", ".fasta"],
                    help="File extensions to scan when --input is a directory.")
    ap.add_argument("--summary_name", default="timing_alignments.tsv",
                    help="Filename for the timing summary TSV in out_dir.")
    args = ap.parse_args()

    clw = which_clustalw()
    if not clw:
        raise SystemExit("Could not find clustalw2 or clustalw in PATH.")

    in_path = Path(args.input)
    out_base = Path(args.out_dir)
    out_base.mkdir(parents=True, exist_ok=True)

# Prepare timing summary
    summary_path = out_base / args.summary_name
    ensure_summary(summary_path)

    # Determine FASTA list
    fasta_list = []
    if in_path.is_dir():
        for ext in args.exts:
            fasta_list += list(in_path.glob(f"*{ext}"))
        if not fasta_list:
            raise SystemExit(f"No FASTA files with {args.exts} in {in_path}")
    else:
        fasta_list = [in_path]

    # Process each FASTA
    for fa in sorted(fasta_list):
        records = list(read_fasta(fa))
        if len(records) == 0:
            print(f"[SKIP] Empty FASTA: {fa}")
            continue
        # Make per-file output dir
        fa_dir = out_base / fa.stem
        fa_dir.mkdir(parents=True, exist_ok=True)

        # Precompute labels for filenames
        labels = [sanitize_label(h) for h, _ in records]

        # Enumerate pairs
        pairs = list(enumerate_pairs(records,
                                     include_self=not args.no_self,
                                     unique=args.unique))

        print(f"[INFO] {fa.name}: {len(records)} sequences → {len(pairs)} alignments")
        for i, j in pairs:
            recA = records[i]
            recB = records[j]
            labA = labels[i]
            labB = labels[j]

            pair_prefix = f"{fa.stem}__{labA}_vs_{labB}"
            pair_fa = fa_dir / f"{pair_prefix}.pair.fa"
            aln_out = fa_dir / f"{pair_prefix}.aln"

            write_pair_fasta(pair_fa, recA, recB)

            t0 = time.monotonic()
            status, note = "OK", ""

            try:
                run_clustalw(clw, pair_fa, aln_out, is_dna=not args.protein, quiet=True)
                elapsed = time.monotonic() - t0
                print(f"[OK] {fa.name}: {labA} vs {labB} -> {aln_out.name}  ({elapsed:.3f}s, {fmt_mmss(elapsed)})")
            except Exception as e:
                elapsed = time.monotonic() - t0
                status, note = "FAIL", str(e).replace("\n", " ")
                print(f"[FAIL] {fa.name}: {labA} vs {labB} -> {note}  ({elapsed:.3f}s, {fmt_mmss(elapsed)})")

            # Record row in summary
            append_summary(summary_path, [
                str(fa),
                str(i),
                str(j),
                labA,
                labB,
                str(pair_fa),
                str(aln_out),
                f"{elapsed:.6f}",
                status,
                note
            ])

    print(f"\nDone. Alignments are in: {out_base.resolve()}")
    print(f"Timing summary: {summary_path.resolve()}")

if __name__ == "__main__":
    main()
