"""The three groups of reference genes that Gan et al. 2025 contrast, taken from published tables only.

    usage:  python analysis/13_gene_groups.py

No gene finding or functional annotation is done here: every gene and every label is read off a
published table.
    coordinates and gene names: the H37Rv RefSeq annotation (NC_000962.3, GCF_000195955.2), which is
        also the coordinate system of the inferred ancestral sequence used for variant calls
    essentiality:  DeJesus et al. 2017 (mBio), Table S3, column "Final Call": ES = essential,
        NE = non-essential. The intermediate calls (ESD, GD, GA, Uncertain) are left out.
    PE/PPE family: membership follows the gene name in the same RefSeq annotation (PE, PPE, PE_PGRS,
        PPE_MPTR and the numbered members), which is how the family is defined in the literature.

The three groups are made disjoint: a PE/PPE gene stays in the PE/PPE group whatever its essentiality
call, so that the comparison between groups is not confounded.

out: results/groups/gene_groups.tsv   one row per reference CDS (locus_tag, coordinates, group)
     results/groups/summary.tsv       genes and coding bases per group
"""
import gzip
import os
import re
import pandas as pd

MTB = os.getcwd()
OUT = f"{MTB}/results/groups"
os.makedirs(OUT, exist_ok=True)

# ---- reference CDS coordinates and gene names
rows = []
with gzip.open(f"{MTB}/data/reference/GCF_000195955.2_ASM19595v2_genomic.gff.gz", "rt") as gff:
    for line in gff:
        if line.startswith("#"):
            continue
        field = line.rstrip("\n").split("\t")
        if len(field) < 9 or field[2] != "CDS":
            continue
        attributes = dict(part.split("=", 1) for part in field[8].split(";") if "=" in part)
        rows.append(dict(locus_tag=attributes.get("locus_tag", ""), start=int(field[3]), end=int(field[4]),
                         strand=field[6], gene=attributes.get("gene", ""),
                         product=attributes.get("product", "").replace("%2C", ",")))
cds = pd.DataFrame(rows).drop_duplicates("locus_tag")
cds["length"] = cds.end - cds.start + 1
print(f"{len(cds)} reference CDS, {cds.length.sum()} coding bases")

# ---- essentiality calls (DeJesus et al. 2017, Table S3)
table_s3 = pd.read_excel(f"{MTB}/data/essential/files/mbo002173137st3.xlsx", sheet_name="ORF Essentiality Calls",
                         skiprows=1)
table_s3.columns = [str(c).strip() for c in table_s3.columns]
orf_column = next(c for c in table_s3.columns if re.search(r"orf|locus|gene.?id", c, re.I))
call_column = next(c for c in table_s3.columns if re.search(r"final.?call", c, re.I))
calls = (table_s3[[orf_column, call_column]].dropna()
         .rename(columns={orf_column: "locus_tag", call_column: "essentiality"}))
calls["locus_tag"] = calls.locus_tag.astype(str).str.strip()
calls["essentiality"] = calls.essentiality.astype(str).str.strip()
calls = calls[calls.locus_tag.str.match(r"^Rv\d")].drop_duplicates("locus_tag")
print(f"essentiality calls: {len(calls)};", calls.essentiality.value_counts().to_dict())
cds = cds.merge(calls, on="locus_tag", how="left")

# ---- groups
is_pe_ppe = cds.gene.str.match(r"^(PE|PPE)(\d|_)", na=False) | cds["product"].str.contains(
    r"\b(PE|PPE)[-_ ]?(family|PGRS|MPTR)", case=False, regex=True, na=False)
cds["group"] = "other"
cds.loc[cds.essentiality == "NE", "group"] = "nonessential"
cds.loc[cds.essentiality == "ES", "group"] = "essential"
cds.loc[is_pe_ppe, "group"] = "pe_ppe"
cds = cds.sort_values("start")
cds["in_all_cds"] = True

# ---- a fifth group: the high-quality parts of the PE/PPE genes, as Gan et al. 2025 define them in
# their Table S3 (segment coordinates in NC_000962.3). Short-read data is unreliable over the
# repetitive remainder of these genes, and they analysed these segments only. Each segment is moved
# inwards to the nearest codon boundary of its own gene so that it stays in frame; these rows are
# extra rows, not part of all_cds, so nothing is counted twice.
table_s3 = pd.read_excel(f"{MTB}/data/gan2025/supp/files/msystems.00898-25-s0010.xlsx", skiprows=1)
table_s3.columns = ["gene", "gene_start", "gene_end", "strand", "segment", "segment_start", "segment_end"]
table_s2 = pd.read_excel(f"{MTB}/data/gan2025/supp/files/msystems.00898-25-s0009.xlsx", skiprows=1)
table_s2.columns = ["seqnames", "start", "end", "width", "strand", "locus_tag", "gene",
                    "include_length", "ratio"]
name_to_locus = dict(zip(table_s2.gene.astype(str).str.strip(), table_s2.locus_tag.astype(str).str.strip()))
coordinates = cds.set_index("locus_tag")[["start", "end", "strand"]].to_dict("index")
segments = []
for _, row in table_s3.iterrows():
    locus_tag = name_to_locus.get(str(row.gene).strip())
    gene = coordinates.get(locus_tag)
    if gene is None or pd.isna(row.segment_start) or pd.isna(row.segment_end):
        continue
    low, high = int(row.segment_start), int(row.segment_end)
    if gene["strand"] == "+":
        offset_low, offset_high = low - gene["start"], high - gene["start"]
    else:
        offset_low, offset_high = gene["end"] - high, gene["end"] - low
    offset_low = max(0, -(-offset_low // 3) * 3)
    length = ((min(offset_high, gene["end"] - gene["start"]) - offset_low + 1) // 3) * 3
    if length < 63:
        continue
    if gene["strand"] == "+":
        start, end = gene["start"] + offset_low, gene["start"] + offset_low + length - 1
    else:
        end, start = gene["end"] - offset_low, gene["end"] - offset_low - length + 1
    segments.append(dict(locus_tag=f"{locus_tag}:{row.segment}", start=start, end=end,
                         strand=gene["strand"], gene=str(row.gene).strip(), product="high-quality segment",
                         length=length, essentiality="", group="pe_ppe_hq", in_all_cds=False))
cds = pd.concat([cds, pd.DataFrame(segments)], ignore_index=True).sort_values(["start", "group"])
print(f"{len(segments)} high-quality PE/PPE segments, {sum(s['length'] for s in segments)} bases, "
      f"from {len({s['gene'] for s in segments})} genes")
cds.to_csv(f"{OUT}/gene_groups.tsv", sep="\t", index=False)

summary = (cds.groupby("group").agg(genes=("locus_tag", "size"), coding_bases=("length", "sum"))
             .reindex(["essential", "nonessential", "pe_ppe", "pe_ppe_hq", "other"]).reset_index())
summary.to_csv(f"{OUT}/summary.tsv", sep="\t", index=False)
print(summary.to_string(index=False))
print("PE/PPE genes with an essentiality call of ES:",
      int(((cds.group == "pe_ppe") & (cds.essentiality == "ES")).sum()))
