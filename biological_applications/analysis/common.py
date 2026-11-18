"""Helpers shared by the analysis scripts: reading a pair's configuration, reading
FASTA / Prodigal / nucmer / GFF files, interval arithmetic, and writing LaTeX number macros.

All scripts are run from the folder that contains README.md. Paths used throughout:
    data/genomes/<accession>/genome.fna, annotation.gff    downloaded genomes (step 1; not in the repository)
    data/genes/<accession>.ffn, .faa                      Prodigal gene calls (step 2; not in the repository)
    work/                                                 intermediate files (not in the repository)
    <pair folder>/results/                                small result tables (in the repository)
"""
import configparser
import os
import pandas as pd

THREADS = os.environ.get("THREADS", "8")      # threads for sourmash, DIAMOND and nucmer


# ---------------------------------------------------------------- configuration of a pair

def read_genomes(pair_dir):
    """genomes.tsv of a pair: one row per genome (role, accession, label, group, key)."""
    return pd.read_csv(f"{pair_dir}/genomes.tsv", sep="\t")


def read_settings(pair_dir):
    """settings.ini of a pair (macro prefix, display names, figure options)."""
    settings = configparser.ConfigParser()
    settings.read(f"{pair_dir}/settings.ini")
    return settings["pair"]


def pair_accessions(pair_dir):
    """Accessions of the two genomes of the pair: A (coordinate reference) and B."""
    genomes = read_genomes(pair_dir)
    a = genomes[genomes.role == "A"].accession.iloc[0]
    b = genomes[genomes.role == "B"].accession.iloc[0]
    return a, b


def work_dir(pair_dir, name):
    """Directory for the intermediate files of one step: work/<pair folder name>/<name>
    (or <pair folder>/<name> if the pair folder itself lies in work/)."""
    pair_dir = os.path.normpath(pair_dir)
    inside_work = pair_dir.split(os.sep)[0] == "work"
    path = os.path.abspath(f"{pair_dir}/{name}" if inside_work else f"work/{os.path.basename(pair_dir)}/{name}")
    os.makedirs(path, exist_ok=True)
    return path


# ---------------------------------------------------------------- sequence files

def read_fasta(path):
    """Return {record name: sequence}. The name is the first word of the header line."""
    sequences = {}
    name = None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            sequences[name] = []
        else:
            sequences[name].append(line)
    return {name: "".join(parts) for name, parts in sequences.items()}


def write_fasta(sequences, path):
    with open(path, "w") as out:
        for name, sequence in sequences.items():
            out.write(f">{name}\n{sequence}\n")


def gc_percent(dna):
    dna = dna.upper()
    return 100 * (dna.count("G") + dna.count("C")) / max(len(dna), 1)


def read_prodigal_genes(path):
    """Gene coordinates from the header lines of a Prodigal FASTA file.
    Prodigal headers look like:  >contig_12 # 3400 # 4600 # -1 # ID=...
    Returns a table with columns gene, contig, start, end, strand."""
    rows = []
    for line in open(path):
        if line.startswith(">"):
            fields = line[1:].split(" # ")
            gene = fields[0].strip()
            contig = gene.rsplit("_", 1)[0]          # the gene name is <contig>_<number>
            rows.append(dict(gene=gene, contig=contig, start=int(fields[1]), end=int(fields[2]),
                             strand="+" if fields[3].strip() == "1" else "-"))
    return pd.DataFrame(rows)


COORDS_COLUMNS = ["start_ref", "end_ref", "start_qry", "end_qry", "len_ref", "len_qry", "identity",
                  "contig_len_ref", "contig_len_qry", "cov_ref", "cov_qry", "contig_ref", "contig_qry"]


def read_coords(path):
    """Alignment blocks written by `show-coords -THrcl` (one row per aligned block)."""
    return pd.read_csv(path, sep="\t", header=None, names=COORDS_COLUMNS)


def read_gff(path):
    """Genes of an NCBI GFF file: contig, start, end, strand, type, locus_tag, product."""
    rows = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] not in ("CDS", "pseudogene", "tRNA", "rRNA", "ncRNA", "tmRNA"):
            continue
        attributes = dict(item.split("=", 1) for item in f[8].split(";") if "=" in item)
        product = attributes.get("product", "pseudogene" if f[2] == "pseudogene" else "")
        rows.append(dict(contig=f[0], start=int(f[3]), end=int(f[4]), strand=f[6], type=f[2],
                         locus_tag=attributes.get("locus_tag", attributes.get("Parent", "")),
                         product=product.replace("%2C", ",").replace("%3B", ";")))
    genes = pd.DataFrame(rows)
    # a CDS split over several lines (or CDS + its pseudogene parent) should be listed once
    return genes.drop_duplicates(["contig", "start", "end"]).reset_index(drop=True)


# ---------------------------------------------------------------- the FracMinHash dN/dS formula

S_OVER_N = 0.77 / 2.23        # ratio of synonymous to non-synonymous sites used in the manuscript


def fmh_dnds_formula(containment_nt, containment_aa, k_aa=7):
    """FracMinHash dN/dS from the nucleotide and the protein containment (k_nt = 3 k_aa)."""
    return S_OVER_N * (1 - containment_aa ** (1 / k_aa)) / (containment_aa ** (1 / k_aa) - containment_nt ** (1 / k_aa))


def fmh_dnds_with_identical_fraction(f, background_nt, background_aa):
    """FracMinHash dN/dS of two genomes of which a fraction f of the coding sequence is identical
    and the remainder has the containments background_nt and background_aa."""
    return fmh_dnds_formula(f + (1 - f) * background_nt, f + (1 - f) * background_aa)


# ---------------------------------------------------------------- intervals

def merge_intervals(intervals):
    """Merge overlapping (start, end) intervals; returns a sorted list of disjoint intervals."""
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def covered_bases(intervals, start, end):
    """Number of positions in [start, end] covered by at least one of the intervals."""
    clipped = [(max(s, start), min(e, end)) for s, e in intervals if e >= start and s <= end]
    return sum(e - s + 1 for s, e in merge_intervals(clipped))


# ---------------------------------------------------------------- LaTeX numbers

def write_numbers(path, prefix, numbers):
    """Write LaTeX macros, one per number quoted in the manuscript.

    numbers is a list of (name, value, meaning). Each becomes
        \\providecommand{\\<prefix><name>}{<value>}   % <meaning>
    LaTeX macro names may contain letters only, hence names like 'FmhDnds' not 'fmh_dnds_1'."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as out:
        out.write("% Generated file -- do not edit by hand. Re-run the analysis instead (see the code repository).\n")
        for name, value, meaning in numbers:
            assert (prefix + name).isalpha(), f"macro name must be letters only: {prefix + name}"
            out.write(f"\\providecommand{{\\{prefix}{name}}}{{{value}}}  % {meaning}\n")
    print(f"wrote {path} ({len(numbers)} numbers)")


def thousands(n):
    """12345 -> '12{,}345' (LaTeX-safe thousands separator)."""
    return f"{int(round(n)):,}".replace(",", "{,}")


def kb(n, digits=0):
    """Base pairs -> kilobases as text, e.g. 38466 -> '38'."""
    return f"{n / 1000:.{digits}f}"
