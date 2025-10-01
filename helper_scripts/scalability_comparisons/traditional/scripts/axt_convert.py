#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Recursively convert pairwise alignments (CLUSTAL .aln or 2-seq FASTA) to AXT,
ensuring the final aligned length is divisible by 3 by appending '-' to BOTH
sequences as needed. Mirrors input directory structure in --out_dir.
Records per-file conversion time to timing_summary.tsv.

NEW: before exiting, concatenates all successfully generated AXT files into
<out_dir>/all_pairs_combined.axt (in sorted, deterministic order).
"""

import argparse
import re
import sys
import time
from pathlib import Path
from typing import List, Tuple, Dict

# ---------- Small utils ----------

def sanitize_label(h: str) -> str:
    """Take first token of header, strip '>', keep safe chars for filenames."""
    token = h[1:].split()[0] if h.startswith('>') else h.split()[0]
    return re.sub(r'[^A-Za-z0-9._|\-]+', '_', token)

def pad_both_to_multiple_of_three(seqA: str, seqB: str) -> Tuple[str, str]:
    """Equalize lengths, then pad both so total length % 3 == 0."""
    if len(seqA) != len(seqB):
        if len(seqA) < len(seqB):
            seqA += "-" * (len(seqB) - len(seqA))
        else:
            seqB += "-" * (len(seqA) - len(seqB))
    r = len(seqA) % 3
    if r != 0:
        pad = 3 - r
        seqA += "-" * pad
        seqB += "-" * pad
    return seqA, seqB

# ---------- FASTA (2-seq) ----------

def read_fasta(path: Path) -> List[Tuple[str, str]]:
    recs: List[Tuple[str, str]] = []
    h, chunks = None, []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith(">"):
                if h is not None:
                    recs.append((h, "".join(chunks)))
                h, chunks = line.strip(), []
            else:
                chunks.append(line.strip())
        if h is not None:
            recs.append((h, "".join(chunks)))
    return recs

# ---------- CLUSTAL (pairwise) ----------

def parse_clustal_pair(path: Path) -> List[Tuple[str, str]]:
    seqs: Dict[str, List[str]] = {}
    order: List[str] = []
    with open(path) as f:
        lines = f.readlines()

    i = 0
    # Skip headers like "CLUSTAL ..." or "MUSCLE ..."
    while i < len(lines) and (not lines[i].strip() or lines[i].upper().startswith("CLUSTAL") or lines[i].startswith("MUSCLE")):
        i += 1

    name_chunk_re = re.compile(r'^(\S+)\s+([A-Za-z\-\.\*]+)')

    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        while i < len(lines) and lines[i].strip():
            m = name_chunk_re.match(lines[i])
            if m:
                name, chunk = m.group(1), m.group(2)
                if name not in seqs:
                    seqs[name] = []
                    order.append(name)
                seqs[name].append(chunk)
            i += 1
        while i < len(lines) and not lines[i].strip():
            i += 1

    if len(order) < 2:
        raise ValueError(f"Expected ≥2 sequences in CLUSTAL file: {path}")
    out = [(nm, "".join(seqs[nm])) for nm in order]
    return out[:2]

# ---------- Core ----------

def load_alignment(path: Path, forced: str = "auto") -> List[Tuple[str, str]]:
    """Return [(label, aligned_seq)] for a pairwise alignment."""
    fmt = forced
    if fmt == "auto":
        if path.suffix.lower() == ".aln":
            fmt = "clustal"
        else:
            with open(path) as f:
                for line in f:
                    if line.strip():
                        fmt = "clustal" if (line.upper().startswith("CLUSTAL") or line.startswith("MUSCLE")) else "fasta"
                        break
    if fmt == "clustal":
        recs = parse_clustal_pair(path)
        return [(sanitize_label(n), s) for n, s in recs]
    elif fmt == "fasta":
        recs = read_fasta(path)
        if len(recs) < 2:
            raise ValueError(f"FASTA must contain at least 2 records: {path}")
        (h1, s1), (h2, s2) = recs[0], recs[1]
        return [(sanitize_label(h1), s1), (sanitize_label(h2), s2)]
    else:
        raise ValueError(f"Unknown format: {forced}")

def write_axt(axt_path: Path, nameA: str, nameB: str, seqA: str, seqB: str) -> None:
    axt_path.parent.mkdir(parents=True, exist_ok=True)
    with open(axt_path, "w") as out:
        out.write(f"{nameA}_vs_{nameB}\n{seqA}\n{seqB}\n\n")

def main():
    ap = argparse.ArgumentParser(
        description="Recursively convert pairwise alignments (.aln CLUSTAL or 2-seq FASTA) to AXT, "
                    "padding to codon length (multiple of 3). Records conversion time and concatenates all AXT."
    )
    ap.add_argument("--input", required=True, help="Alignment file OR directory (recursively scanned).")
    ap.add_argument("--out_dir", required=True, help="Where to write .axt, timing_summary.tsv, and all_pairs_combined.axt")
    ap.add_argument("--format", choices=["auto","clustal","fasta"], default="auto", help="Force input format (default: auto)")
    ap.add_argument("--exts", nargs="+", default=[".aln",".fa",".fna",".fasta"], help="Extensions to scan recursively if --input is a directory")
    args = ap.parse_args()

    in_path = Path(args.input)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Build list of inputs (recurse if directory)
    if in_path.is_dir():
        inputs: List[Path] = []
        for ext in args.exts:
            inputs.extend(sorted(in_path.rglob(f"*{ext}")))
        if not inputs:
            sys.exit(f"No input files with {args.exts} under {in_path}")
    else:
        inputs = [in_path]

    summary = out_dir / "timing_axt.tsv"
    with open(summary, "w") as sf:
        sf.write("input_file\taxt_file\tconvert_seconds\tstatus\tnote\n")

    total, n_ok = 0.0, 0
    generated_axts: List[Path] = []

    for aln in inputs:
        status, note = "OK", ""
        t0 = time.monotonic()
        try:
            recs = load_alignment(aln, forced=args.format)
            (nameA, sA), (nameB, sB) = recs[0], recs[1]
            sA2, sB2 = pad_both_to_multiple_of_three(sA, sB)

            if in_path.is_dir():
                rel = aln.relative_to(in_path)
                axt_path = (out_dir / rel).with_suffix(".axt")
            else:
                axt_path = (out_dir / aln.name).with_suffix(".axt")

            write_axt(axt_path, nameA, nameB, sA2, sB2)
            generated_axts.append(axt_path)

            dt = time.monotonic() - t0
            total += dt; n_ok += 1
            print(f"[OK] {aln} -> {axt_path} ({dt:.3f}s)")
            with open(summary, "a") as sf:
                sf.write(f"{aln}\t{axt_path}\t{dt:.6f}\t{status}\t{note}\n")
        except Exception as e:
            dt = time.monotonic() - t0
            status, note = "FAIL", str(e).replace("\n"," ")
            print(f"[FAIL] {aln}: {note}")
            with open(summary, "a") as sf:
                sf.write(f"{aln}\t\t\t{status}\t{note}\n")

    # ---- Concatenate all generated AXT files (sorted) ----
    combined_path = out_dir / "sequences.axt"
    if generated_axts:
        generated_axts = sorted(set(p.resolve() for p in generated_axts))
        with open(combined_path, "w") as out:
            for p in generated_axts:
                with open(p, "r") as inp:
                    out.write(inp.read())
        print(f"[OK] Concatenated {len(generated_axts)} AXT files -> {combined_path}")
    else:
        print("[WARN] No AXT files generated; skipping concatenation.")

    print("\n==== Conversion Timing ====")
    avg = (total / n_ok) if n_ok else 0.0
    print(f"Converted: {n_ok} files, total {total:.3f}s, avg {avg:.3f}s")
    print(f"Summary: {summary.resolve()}")

if __name__ == "__main__":
    main()
