#!/usr/bin/env python
"""Step 18. Supplementary data files of the manuscript.

    usage:  python analysis/18_supplementary_data.py

Output (in output/data/):
    S1_Data_genome_pairs.tsv   the 140 genome pairs of the hierarchical edge bundling figure with
                               their FracMinHash dN/dS, the location of their shared k-mers and
                               their class (all_pairs/results/all_pairs_summary.tsv)
    S2_Data_genomes.tsv        every genome assembly used: the 243 genomes of the 140 pairs and
                               the additional genomes compared with the pairs examined in detail
"""
import os
import pandas as pd

os.makedirs("output/data", exist_ok=True)
pd.read_csv("all_pairs/results/all_pairs_summary.tsv", sep="\t").to_csv("output/data/S1_Data_genome_pairs.tsv", sep="\t", index=False)

rows = []
gtdb = pd.read_csv("all_pairs/genomes.tsv", sep="\t")
for genome in gtdb.itertuples():
    rows.append(dict(accession=genome.accession, used_in="hierarchical edge bundling figure", name=genome.ncbi_organism_name,
                     gtdb_species=genome.gtdb_species, genome_type=genome.genome_type))
for pair in pd.read_csv("pairs/examined_pairs.tsv", sep="\t").itertuples():
    for genome in pd.read_csv(f"{pair.folder}/genomes.tsv", sep="\t").itertuples():
        role = {"A": "genome of the pair", "B": "genome of the pair"}.get(genome.role, f"compared genome ({genome.group})")
        rows.append(dict(accession=genome.accession, used_in=f"pair {pair.number} ({pair.name}): {role}", name=genome.label))
table = pd.DataFrame(rows)
# one row per accession; a genome that is used in several places lists all of them
merged = table.groupby("accession", sort=False).agg(dict(used_in="; ".join, name="first", gtdb_species="first", genome_type="first")).reset_index()
merged.to_csv("output/data/S2_Data_genomes.tsv", sep="\t", index=False)
print(f"{len(merged)} genome assemblies")
