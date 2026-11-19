"""The isolates of the equivalence check (step 31), as input lists of the FracMinHash dN/dS tool.

    usage:  python checks/31a_check_inputs.py [out_folder]

1,000 assembly-derived isolates chosen at random (seed 20261118), plus the reference, for the PE/PPE
high-quality segments and the whole coding genome; and the first 80 of them for the whole coding
genome, because main mode cannot read back its own sketch file of more whole coding genomes (step 31).

out (work/equivalence/ by default): pe_ppe_hq_1000.csv, all_cds_1000.csv, all_cds_80.csv
"""
import os
import random
import sys

MTB = os.getcwd()
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{MTB}/work/equivalence"
os.makedirs(OUT, exist_ok=True)
names = [line.strip() for line in open(f"{MTB}/work/fasta/assembly_names.txt") if line.strip()]
random.seed(20261118)
chosen = sorted(random.sample(names, 1000))


def write(group, isolates, path):
    rows = ["name,genome_filename,protein_filename"] + [
        f"{n},{MTB}/work/fasta/{group}/{n}.ffn,{MTB}/work/fasta/{group}/{n}.faa"
        for n in ["ancestral_reference"] + isolates]
    open(path, "w").write("\n".join(rows) + "\n")


for group in ("all_cds", "pe_ppe_hq"):
    write(group, chosen, f"{OUT}/{group}_1000.csv")
write("all_cds", chosen[:80], f"{OUT}/all_cds_80.csv")
print(f"{len(chosen)} isolates; lists in {OUT}")
