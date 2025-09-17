#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import random
import argparse
from pathlib import Path

# -----------------------
# CONFIG (edit as needed)
# -----------------------
wd = '/data/jzr5814/sourmash_dnds_estimation/tests/data/ecoli_10strains_refseq_redownload'

strain1 = f'{wd}/GCF_027925825.1/cds_from_genomic.fna'
strain2 = f'{wd}/GCF_027925805.1/cds_from_genomic.fna'
strain3 = f'{wd}/GCF_024300685.1/cds_from_genomic.fna'

sample = 1000                 # how many gene sets to sample
min_len = 1000             # length filter applied to strain1 only
seed = 42                  # set None for non-deterministic sampling

# Base output dir
base_out = Path("/data/jzr5814/dnds_scalability_comparison")

# -----------------------
# Helpers
# -----------------------

# All runtime parameters now come from CLI args (see main()).

CORE_ID_RE = re.compile(r'(cds_WP_\d+\.\d)')

def read_fasta(path):
    """Yield (header, seq) from FASTA file."""
    header, seq_chunks = None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_chunks)
                header = line.strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(seq_chunks)

def wrap_seq(s, width=60):
    return '\n'.join(s[i:i+width] for i in range(0, len(s), width))

def core_id_from_header(header):
    """Extract cds_WP_... core id from a FASTA header; return None if absent."""
    m = CORE_ID_RE.search(header)
    return m.group(1) if m else None

def index_by_core_id(fasta_path):
    """
    Build dict core_id -> (header, seq) for a .cds.fna file.
    If multiple entries share a core_id, last one wins.
    """
    idx = {}
    for h, s in read_fasta(fasta_path):
        cid = core_id_from_header(h)
        if cid:
            idx[cid] = (h, s)
    return idx

# -----------------------
# Main workflow
# -----------------------

def main():
    # Make subfolder named "sample_<N>"
    out_dir = base_out / f"sample_{sample}"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Index all three strains by core id
    idx1 = index_by_core_id(strain1)
    idx2 = index_by_core_id(strain2)
    idx3 = index_by_core_id(strain3)

    # Apply length filter to strain1 only
    eligible_in_1 = {cid for cid, (h, s) in idx1.items() if len(s) > min_len}

    # Keep only core IDs present in ALL THREE strains
    common = eligible_in_1 & idx2.keys() & idx3.keys()

    if not common:
        print(f"No common CDS (length>{min_len}) across all three strains. Nothing to do.")
        return

    # Sample
    ids = sorted(common)
    if seed is not None:
        random.seed(seed)
    if len(ids) < sample:
        print(f"Only {len(ids)} common CDS available; reducing sample from {sample} to {len(ids)}.")
        sample_n = len(ids)
    else:
        sample_n = sample

    picked = random.sample(ids, sample_n)

    # Write one FASTA per picked core id, containing the three sequences
    manifest_lines = ["core_id\toutfile\tstrain_file\theader"]
    for cid in picked:
        out_fa = out_dir / f"{cid}.fasta"
        with open(out_fa, "w") as out:
            for src_path, idx in [(strain1, idx1), (strain2, idx2), (strain3, idx3)]:
                header, seq = idx[cid]
                h = header if header.startswith('>') else f">{header}"
                out.write(f"{h}\n{wrap_seq(seq)}\n")
                manifest_lines.append(f"{cid}\t{out_fa.name}\t{Path(src_path).name}\t{header}")

        print(f"Wrote {out_fa}")

    # Also write a manifest TSV for bookkeeping
    manifest_path = out_dir / "manifest.tsv"
    with open(manifest_path, "w") as mf:
        mf.write("\n".join(manifest_lines) + "\n")
    print(f"Wrote manifest: {manifest_path}")

if __name__ == "__main__":
    main()
