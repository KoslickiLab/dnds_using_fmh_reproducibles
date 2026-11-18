#!/usr/bin/env python
"""Re-create the two input tables of the 140-pair analysis from their public sources.

    usage:  python all_pairs/make_inputs.py          (downloads about 120 MB; the tables ship with the repository)

all_pairs/pairs.tsv
    The genome pairs of the manuscript's hierarchical edge bundling figure, with their FracMinHash
    dN/dS: figures/connect_species.csv of KoslickiLab/DnDs-visualization (commit 060896e), one row
    per pair, in the same order; a pair identifier P001, P002, ... is added.
all_pairs/genomes.tsv
    GTDB release 214 (214.1) metadata of the genomes in these pairs: species and taxonomy, genome
    type (isolate, MAG or SAG, from the NCBI genome category), CheckM completeness and contamination,
    number of contigs and genome size.

The downloads are written to work/inputs/.
"""
import os
import subprocess
import pandas as pd

CONNECT = ("https://raw.githubusercontent.com/KoslickiLab/DnDs-visualization/"
           "060896ed9ba86eddb1278b9bfae4c677cd044ec8/figures/connect_species.csv")
GTDB = "https://data.gtdb.ecogenomic.org/releases/release214/214.1"
GENOME_TYPE = {"none": "isolate", "derived from metagenome": "MAG",
               "derived from environmental sample": "MAG", "derived from single cell": "SAG"}

os.makedirs("work/inputs", exist_ok=True)


def download(url, path):
    if not os.path.exists(path):
        subprocess.run(["curl", "-sSfL", "-o", path, url], check=True)
    return path


# ---- the pairs
connect = pd.read_csv(download(CONNECT, "work/inputs/connect_species.csv"))
pairs = pd.DataFrame(dict(pair=[f"P{i:03d}" for i in range(1, len(connect) + 1)], genome_1=connect["from"],
                          genome_2=connect["to"], fmh_dnds_in_paper=connect.dndsvalue))
pairs.to_csv("all_pairs/pairs.tsv", sep="\t", index=False)

# ---- GTDB R214 metadata of their genomes (GTDB writes accessions as RS_GCF_... or GB_GCA_...)
wanted = set(pairs.genome_1) | set(pairs.genome_2)
tables = []
for domain in ("ar53", "bac120"):
    path = download(f"{GTDB}/{domain}_metadata_r214.tsv.gz", f"work/inputs/{domain}_metadata_r214.tsv.gz")
    metadata = pd.read_csv(path, sep="\t", low_memory=False)
    metadata["accession"] = metadata.accession.str[3:]
    tables.append(metadata[metadata.accession.isin(wanted)])
metadata = pd.concat(tables)
genomes = pd.DataFrame(dict(
    accession=metadata.accession,
    gtdb_species=metadata.gtdb_taxonomy.str.split(";s__").str[1],
    gtdb_taxonomy=metadata.gtdb_taxonomy,
    genome_type=metadata.ncbi_genome_category.map(GENOME_TYPE),
    ncbi_organism_name=metadata.ncbi_organism_name,
    checkm_completeness=metadata.checkm_completeness,
    checkm_contamination=metadata.checkm_contamination,
    contig_count=metadata.contig_count,
    genome_size=metadata.genome_size)).sort_values("accession")
assert len(genomes) == len(wanted), "some genomes of the pairs are missing from the GTDB R214 metadata"
genomes.to_csv("all_pairs/genomes.tsv", sep="\t", index=False)
print(f"wrote all_pairs/pairs.tsv ({len(pairs)} pairs) and all_pairs/genomes.tsv ({len(genomes)} genomes)")
