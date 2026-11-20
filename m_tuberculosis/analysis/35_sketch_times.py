"""How long the tool took to sketch each run, from its own log lines (logs/32_run_reference_mode.log).

    usage:  python analysis/35_sketch_times.py

The tool logs "Sketching data fasta file: <shard csv>" when it starts a shard and "Successfully
sketched data in <s> seconds" when the shard is done. A run's sketching time is from its first start
to its last finish; the shards of one run are sketched side by side, one core each. Shards and genomes are counted from the shard folders.

out: results/tool/sketching.tsv   run, shards, genomes, start, end, wall seconds, core-hours
"""
import glob
import os
import re
from datetime import datetime

import pandas as pd

MTB = os.getcwd()
STAMP = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d\.\d+)")
START = re.compile(r"Sketching data fasta file: (\S+/work/tool/([^/]+)/sketch_shards/\d+/genomes\.csv)")
FINISH = re.compile(r"Successfully sketched data in ([0-9.]+) seconds")

runs = {}
current = None
for line in open(f"{MTB}/logs/32_run_reference_mode.log"):
    stamp = STAMP.match(line)
    if not stamp:
        continue
    when = datetime.strptime(stamp.group(1), "%Y-%m-%d %H:%M:%S.%f")
    started = START.search(line)
    if started:
        csv, current = started.groups()
        run = runs.setdefault(current, dict(start=when, end=when, shards=set()))
        run["start"] = min(run["start"], when)
        run["shards"].add(csv)
    elif FINISH.search(line) and current:
        runs[current]["end"] = max(runs[current]["end"], when)

rows = []
for name, run in runs.items():
    # shards and genomes from the shard folders (with hundreds of shards logging at once, a few log
    # lines interleave and do not parse; the first start and last finish are not affected)
    shard_csvs = glob.glob(f"{MTB}/work/tool/{name}/sketch_shards/*/genomes.csv")
    genomes = sum(sum(1 for _ in open(csv)) - 1 for csv in shard_csvs)
    seconds = (run["end"] - run["start"]).total_seconds()
    rows.append(dict(run=name, shards=len(shard_csvs), genomes=genomes, start=run["start"], end=run["end"],
                     wall_seconds=round(seconds, 1), core_hours=round(seconds * len(shard_csvs) / 3600, 1)))
table = pd.DataFrame(rows)
os.makedirs(f"{MTB}/results/tool", exist_ok=True)
table.to_csv(f"{MTB}/results/tool/sketching.tsv", sep="\t", index=False)
print(table.to_string(index=False))
