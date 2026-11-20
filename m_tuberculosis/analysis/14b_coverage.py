"""How much of the reference each assembly covers (from step 14's alignments).

    usage:  python analysis/14b_coverage.py

out: results/inputs/coverage.tsv   per assembly: bases of the ancestral sequence (4,411,532 bp) that its
                                   alignment covers; positions not covered keep the reference base
"""
import glob
import os

import pandas as pd

MTB = os.getcwd()
rows = [dict(sample=os.path.basename(path)[:-len(".cov")], covered_bases=int(open(path).read().strip() or 0))
        for path in sorted(glob.glob(f"{MTB}/work/variants/*.cov"))]
table = pd.DataFrame(rows)
table.to_csv(f"{MTB}/results/inputs/coverage.tsv", sep="\t", index=False)
print(f"{len(table)} assemblies; mean coverage {100 * table.covered_bases.mean() / 4411532:.1f}% of the reference")
