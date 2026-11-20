#!/usr/bin/env python
"""Step 8. Scaffolds are contigs joined across gaps (runs of N). A shared region that is
separated from the rest of its scaffold by gaps on both sides is NOT proven to be integrated
there. This script lists, for each shared region, the gaps inside and within 20 kb of it, in
genome A and in genome B.

    usage:  python analysis/08_scaffold_gaps.py <pair folder>
    output: <pair folder>/results/gaps_near_regions.tsv
"""
import sys
import re
import pandas as pd
from common import pair_accessions, read_fasta

WINDOW = 20000
pair_dir = sys.argv[1]
accession_a, accession_b = pair_accessions(pair_dir)
blocks = pd.read_csv(f"{pair_dir}/results/blocks.tsv", sep="\t")
regions = pd.read_csv(f"{pair_dir}/results/regions.tsv", sep="\t")


def gaps(sequence):
    """(start, end) of every run of at least 10 N, 1-based."""
    return [(m.start() + 1, m.end()) for m in re.finditer("N{10,}", sequence.upper())]


rows = []
genome = {"A": read_fasta(f"data/genomes/{accession_a}/genome.fna"), "B": read_fasta(f"data/genomes/{accession_b}/genome.fna")}
gap_table = {side: {contig: gaps(seq) for contig, seq in genome[side].items()} for side in "AB"}
for region in regions.itertuples():
    mine = blocks[blocks.region == region.region]
    # the stretch of each genome-B contig that this region aligns to
    stretches = [("A", region.contig, region.start, region.end)]
    for contig, part in mine.groupby("contig_qry"):
        stretches.append(("B", contig, int(part[["start_qry", "end_qry"]].min().min()), int(part[["start_qry", "end_qry"]].max().max())))
    for side, contig, start, end in stretches:
        contig_length = len(genome[side][contig])
        left = [g for g in gap_table[side][contig] if start - WINDOW <= g[1] < start]
        inside = [g for g in gap_table[side][contig] if g[0] >= start and g[1] <= end]
        right = [g for g in gap_table[side][contig] if end < g[0] <= end + WINDOW]
        rows.append(dict(region=region.region, genome=side, contig=contig, contig_length=contig_length, start=start, end=end,
                         bp_to_contig_start=start - 1, bp_to_contig_end=contig_length - end,
                         gaps_left_20kb=len(left), gaps_inside=len(inside), gaps_right_20kb=len(right),
                         nearest_gap_left_bp=start - max([g[1] for g in left], default=float("nan")),
                         nearest_gap_right_bp=min([g[0] for g in right], default=float("nan")) - end))
table = pd.DataFrame(rows)
table.to_csv(f"{pair_dir}/results/gaps_near_regions.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(table.to_string(index=False))
