#!/usr/bin/env python
"""Step 9. Is each shared region present in other genomes?

    usage:  python analysis/09_presence_in_other_genomes.py <pair folder>

Control genomes (rows of genomes.tsv with role "control") answer two different questions:
  * independent assemblies of the SAME strain: if a region is there too, it is part of that
    strain's genome and not a contaminant of one sequencing run;
  * OTHER strains and species: if a region is missing while the DNA on either side of it is
    present, the region is an insertion (mobile), not ancestral sequence.

For every genome other than A, and every region, we report
    covered        fraction of the region covered by alignment blocks of >= 90% identity
    identity       length-weighted mean identity of those blocks
    flank_covered  fraction of the 10 kb on either side of the region that is aligned at all
    located_on     the contig of the other genome that holds most of the region, and its length

Output: <pair folder>/results/presence.tsv
"""
import sys
import pandas as pd
from common import read_genomes, read_coords, covered_bases, work_dir

MIN_IDENTITY = 90.0    # percent, for calling a region present
FLANK = 10000          # bp on each side

pair_dir = sys.argv[1]
results = f"{pair_dir}/results"
genomes = read_genomes(pair_dir)
others = genomes[genomes.role != "A"]
regions = pd.read_csv(f"{results}/regions.tsv", sep="\t")

rows = []
for genome in others.itertuples():
    blocks = read_coords(f"{work_dir(pair_dir, 'alignments')}/{genome.accession}.coords")
    for region in regions.itertuples():
        on_contig = blocks[blocks.contig_ref == region.contig]
        good = on_contig[on_contig.identity >= MIN_IDENTITY]
        inside = good[(good.end_ref >= region.start) & (good.start_ref <= region.end)]
        covered = covered_bases(list(zip(good.start_ref, good.end_ref)), region.start, region.end)
        everything = list(zip(on_contig.start_ref, on_contig.end_ref))
        left = covered_bases(everything, region.start - FLANK, region.start - 1)
        right = covered_bases(everything, region.end + 1, region.end + FLANK)
        identity = (inside.identity * inside.len_ref).sum() / inside.len_ref.sum() if len(inside) else float("nan")
        # where the region sits in the other genome: the contig that holds most of it, and that contig's length
        if len(inside):
            aligned_per_contig = inside.groupby("contig_qry").len_ref.sum()
            located_on = aligned_per_contig.idxmax()
            located_on_length = int(inside[inside.contig_qry == located_on].contig_len_qry.iloc[0])
        else:
            located_on, located_on_length = "", 0
        rows.append(dict(accession=genome.accession, region=region.region,
                         covered=round(covered / region.span_bp, 4), identity=round(identity, 2),
                         flank_covered=round((left + right) / (2 * FLANK), 4),
                         located_on=located_on, located_on_length=located_on_length))
table = pd.DataFrame(rows)
table.to_csv(f"{results}/presence.tsv", sep="\t", index=False)

# a compact view: one row per genome, one column per region, values = fraction of the region covered
wide = table.pivot(index="accession", columns="region", values="covered")
wide = wide[sorted(wide.columns, key=lambda r: int(r[1:]))].loc[others.accession]
pd.set_option("display.width", 250)
print(wide.round(2).to_string())
