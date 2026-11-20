"""The ratio of synonymous to non-synonymous sites in this organism's coding sequences.

    usage:  python analysis/24_site_ratio.py

The estimator multiplies PdN/PdS by the constant 0.77/2.23, which is the number of synonymous and
non-synonymous sites per codon under uniform codon usage. M. tuberculosis is 65% GC and its codon
usage is not uniform, so the constant that belongs to these sequences is not exactly that one. This
counts the sites the way the counting methods do (for each position of each codon, the share of the
three possible substitutions that do not change the amino acid) and reports what the constant would
be, per gene group.

out: results/site_ratio.tsv
"""
import os
import sys
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import anchored

MTB = os.getcwd()
CODONS = anchored.CODONS
GROUPS = anchored.GROUPS

synonymous_sites = {}
for codon, residue in CODONS.items():
    count = 0.0
    for position in range(3):
        for base in "ACGT":
            if base == codon[position]:
                continue
            other = codon[:position] + base + codon[position + 1:]
            if CODONS[other] == residue:
                count += 1 / 3
    synonymous_sites[codon] = count

backbone = anchored.read_backbone()
genes = anchored.gene_table()
rows = []
totals = {group: [0.0, 0.0] for group in GROUPS}
for locus_tag, group, in_all, cds, protein in anchored.coding_sequences(backbone, genes):
    synonymous = non_synonymous = 0.0
    for i in range(0, len(cds) - len(cds) % 3, 3):
        codon = cds[i:i + 3]
        if codon in synonymous_sites:
            synonymous += synonymous_sites[codon]
            non_synonymous += 3 - synonymous_sites[codon]
    for target in (["all_cds"] if in_all else []) + ([group] if group in GROUPS else []):
        totals[target][0] += synonymous
        totals[target][1] += non_synonymous

for group, (synonymous, non_synonymous) in totals.items():
    if synonymous == 0:
        continue
    codons = (synonymous + non_synonymous) / 3
    rows.append(dict(group=group, codons=int(codons),
                     synonymous_sites_per_codon=synonymous / codons,
                     non_synonymous_sites_per_codon=non_synonymous / codons,
                     constant_for_these_sequences=synonymous / non_synonymous,
                     constant_used=0.77 / 2.23,
                     ratio_to_constant_used=(synonymous / non_synonymous) / (0.77 / 2.23)))
table = pd.DataFrame(rows)
table.to_csv(f"{MTB}/results/site_ratio.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(table.round(4).to_string(index=False))
