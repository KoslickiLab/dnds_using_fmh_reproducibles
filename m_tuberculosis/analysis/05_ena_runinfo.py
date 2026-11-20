"""ENA run metadata (layout, platform, bytes, sample) for all 49,954 runs of Gan et al. Table S1.

    usage:  python analysis/05_ena_runinfo.py
    out:  data/isolates/ena_runinfo.tsv
"""
import pandas as pd, requests, io, time
runs = pd.read_csv("data/isolates/gan_runs.tsv", sep="\t").run.drop_duplicates().tolist()
fields = "run_accession,sample_accession,library_layout,library_strategy,library_source,instrument_platform,instrument_model,read_count,base_count,fastq_bytes,sra_bytes,tax_id,first_public"
out = []
for i in range(0, len(runs), 400):
    chunk = runs[i:i + 400]
    for attempt in range(5):
        r = requests.post("https://www.ebi.ac.uk/ena/portal/api/search",
                          data={"result": "read_run", "includeAccessions": ",".join(chunk), "fields": fields, "format": "tsv", "limit": 0}, timeout=300)
        if r.ok:
            break
        time.sleep(10)
    r.raise_for_status()
    out.append(pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str))
df = pd.concat(out)
df.to_csv("data/isolates/ena_runinfo.tsv", sep="\t", index=False)
print(len(runs), "runs requested;", len(df), "returned")
print(df.library_layout.value_counts().to_string()); print(df.instrument_platform.value_counts().to_string())
