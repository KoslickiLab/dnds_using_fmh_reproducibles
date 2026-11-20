#!/usr/bin/env python
"""Step 17 (optional; needs step 16). Read depth and read pairs at the region shared by the
archaeal and the bacterial metagenome-assembled genome.

    usage:  python analysis/17_read_level_summary.py pairs/methanobrevibacter_saccharibacteria

If the shared region belonged to both organisms, the reads of both would pile up on its single
copy in the reference, and its depth would be the sum of the depths of the two genomes. If it
belongs to the bacterium only, its depth equals that of the bacterial sequence around it.
Likewise, read pairs that start at an end of the region continue into the DNA that lies next
to it in the sample: bacterial flanks, or archaeal contigs. Only reads within END_WINDOW of an
end can belong to such pairs (the inserts of this library are shorter than 500 bp). Pairs with
a read in the interior of the region and the other read elsewhere are chimeric pairs, which
every library contains; they show how many such pairs to expect by chance near the ends.

Output (in <pair folder>/results/):
    read_depth.tsv         median read depth of the shared region, of the rest of its scaffold
                           (the bacterial flanks), of the bacterial and of the archaeal genome,
                           and the depth expected if both organisms carried the region
    read_depth_windows.tsv median depth in 1-kb windows along the scaffold of genome A that
                           carries the region (for the figure)
    read_pairs.tsv         read pairs with one read inside the region and the other outside it,
                           by the position of the inside read (at an end / interior) and by
                           where the other read lies
"""
import sys
import subprocess
import io
import pandas as pd
from common import pair_accessions, read_fasta, work_dir

MIN_MAPQ = 2       # the read outside the region must be placed unambiguously
END_WINDOW = 500   # bp; reads this close to an end of the region can belong to a pair that spans the end
READ_LENGTH = 150  # bp
WINDOW = 1000      # bp

pair_dir = sys.argv[1]
results = f"{pair_dir}/results"
bam = f"{work_dir(pair_dir, 'reads')}/mapped.bam"
accession_a, accession_b = pair_accessions(pair_dir)
contigs_a = set(read_fasta(f"data/genomes/{accession_a}/genome.fna"))
region = pd.read_csv(f"{results}/regions.tsv", sep="\t").iloc[0]        # the pair has one shared region

# ---- 1. depth at every position of the reference (samtools depth -a: positions without reads included)
depth = pd.read_csv(io.StringIO(subprocess.run(["samtools", "depth", "-a", bam], capture_output=True, text=True, check=True).stdout),
                    sep="\t", header=None, names=["contig", "position", "depth"])
on_scaffold = depth[depth.contig == region.contig]
in_region = on_scaffold[(on_scaffold.position >= region.start) & (on_scaffold.position <= region.end)]
flanks = on_scaffold[(on_scaffold.position < region.start) | (on_scaffold.position > region.end)]
genome_a = depth[depth.contig.isin(contigs_a) & ~depth.index.isin(in_region.index)]
genome_b = depth[~depth.contig.isin(contigs_a)]
rows = [("shared region", in_region.depth.median()),
        ("rest of the scaffold that carries the region (genome A)", flanks.depth.median()),
        ("genome A without the region", genome_a.depth.median()),
        ("genome B (contigs that consist of shared sequence left out)", genome_b.depth.median()),
        ("expected for the region if both genomes carried it", flanks.depth.median() + genome_b.depth.median())]
table = pd.DataFrame(rows, columns=["sequence", "median_depth"])
table.to_csv(f"{results}/read_depth.tsv", sep="\t", index=False)
print(table.to_string(index=False))

windows = on_scaffold.assign(window_start=(on_scaffold.position - 1) // WINDOW * WINDOW + 1)
windows = windows.groupby("window_start").depth.median().reset_index().rename(columns={"depth": "median_depth"})
windows.insert(0, "contig", region.contig)
windows.to_csv(f"{results}/read_depth_windows.tsv", sep="\t", index=False)

# ---- 2. read pairs with one read inside the region and the other read elsewhere
columns = ["name", "flag", "contig", "position", "mapq", "cigar", "mate_contig", "mate_position"]
SKIP = "0x90C"     # skip unmapped reads, reads with an unmapped mate, secondary and supplementary alignments


def alignments(*arguments):
    """Rows of `samtools view` as a table (first eight SAM columns)."""
    text = subprocess.run(["samtools", "view", "-F", SKIP, *arguments], capture_output=True, text=True, check=True).stdout
    return pd.read_csv(io.StringIO(text), sep="\t", header=None, usecols=range(8), names=columns)


def inside_region(contig, position):
    return (contig == region.contig) & (position >= region.start - READ_LENGTH) & (position <= region.end)


# reads that overlap the region and whose mate was mapped outside it
reads = alignments(bam, f"{region.contig}:{region.start}-{region.end}")
mate_contig = reads.mate_contig.where(reads.mate_contig != "=", reads.contig)
names = reads.name[~inside_region(mate_contig, reads.mate_position)].drop_duplicates()
names.to_csv(f"{work_dir(pair_dir, 'reads')}/pairs_leaving_region.txt", index=False, header=False)
# both reads of those pairs: the read inside the region and the "other" read outside it
both = alignments("-N", f"{work_dir(pair_dir, 'reads')}/pairs_leaving_region.txt", bam)
is_inside = inside_region(both.contig, both.position)
pairs = both[is_inside].merge(both[~is_inside], on="name", suffixes=("_inside", "_other"))
pairs = pairs[pairs.mapq_other >= MIN_MAPQ]
distance_to_end = pd.concat([pairs.position_inside - region.start, region.end - (pairs.position_inside + READ_LENGTH)], axis=1).min(axis=1)
pairs["inside_read"] = ["at an end of the region" if d <= END_WINDOW else "interior of the region" for d in distance_to_end]


def where(contig):
    if contig == region.contig:
        return "flanking sequence of the same scaffold (genome A)"
    return "another contig of genome A" if contig in contigs_a else "a contig of genome B"


pairs["other_read_lies_on"] = pairs.contig_other.map(where)
places = ["flanking sequence of the same scaffold (genome A)", "another contig of genome A", "a contig of genome B"]
counts = pairs.groupby(["inside_read", "other_read_lies_on"]).size()
table = pd.DataFrame([dict(inside_read=position, other_read_lies_on=place, read_pairs=int(counts.get((position, place), 0)))
                      for position in ("at an end of the region", "interior of the region") for place in places])
table.to_csv(f"{results}/read_pairs.tsv", sep="\t", index=False)
print(table.to_string(index=False))
