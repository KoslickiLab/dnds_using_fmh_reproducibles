# Biological applications of FracMinHash dN/dS

Code, configuration and result tables for the section "Biological application of FracMinHash
dN/dS" of the manuscript *Leveraging FracMinHash Containment for Genomic dN/dS*.

A FracMinHash dN/dS above 1 between two distantly related genomes indicates that the genomes
share near-identical DNA. This repository locates that DNA for the 140 genome pairs of the
manuscript's hierarchical edge bundling figure and examines seven pairs in detail: one pair of
metagenome-assembled genomes (a contig assigned to the wrong genome) and six pairs of isolate
genomes (mobile genetic elements shared by distantly related bacteria).

The method itself is in [KoslickiLab/dnds-using-fmh](https://github.com/KoslickiLab/dnds-using-fmh)
(used here unmodified, at commit `ee94eb7`). The other analyses of the manuscript (simulations,
scalability, comparison with published values) are described in the [README of this repository](../README.md).

## Reproduce everything

```bash
conda env create -f environment.yml
conda activate fmh_dnds_applications
bash reproduce_all.sh
```

`reproduce_all.sh` downloads about 280 genomes from NCBI (1 GB), runs every analysis step, draws
every figure and finishes with a check that the gene calls and the FracMinHash dN/dS values equal
those behind the manuscript. By default it analyses 4 genome pairs at a time with 2 threads each;
`JOBS=40 bash reproduce_all.sh` runs 40 at a time. On our Linux server the complete run took
20 minutes with `JOBS=40`, including the downloads; the largest genome pair needs 1.5 minutes and
1.8 GB of memory, so the defaults should finish in about an hour. All scripts are run from the
folder that contains this README.

The results of our run are part of the repository (see "What is in the repository"). A complete
run reproduces every table, every number file and every figure byte for byte, so the regenerated
files can be compared with ours with `git status` or `git diff`.

Only the figures:

```bash
bash figures/make_all_figures.sh      # one minute; reads the result tables, needs no genomes
```

## Figures

Each figure of the manuscript is drawn by one script from the result tables in this repository.

| Figure in the manuscript | Script | Output in `output/figures/` |
|---|---|---|
| A bacterial contig assigned to an archaeal MAG | `figures/figure_mag_misbinning.py` | `bioapp_mag.pdf` |
| *F. duncaniae* A2-165 and *R. inulinivorans* A2-194 | `figures/figure_shared_regions.py pairs/faecalibacterium_roseburia` | `bioapp_faecalibacterium_roseburia.pdf` |
| *C. massiliensis* and *B. intestinihominis* | `figures/figure_shared_regions.py pairs/bacteroidales` | `bioapp_bacteroidales.pdf` |
| *M. vaginalis* and *A. vaginalis* | `figures/figure_shared_regions.py pairs/peptoniphilaceae_1` | `bioapp_peptoniphilaceae.pdf` |
| *R. intestinalis* L1-82 and *R. zhanii* BX1005 | `figures/figure_shared_regions.py pairs/roseburia` | `bioapp_roseburia.pdf` |
| *L. pacificus* and *S. pentaromativorans* | `figures/figure_shared_regions.py pairs/rhodobacteraceae` | `bioapp_rhodobacteraceae.pdf` |
| Supplement: near-identical genes in the 140 genome pairs (box plots) | `figures/figure_supp_near_identical_genes.py` | `bioapp_supp_near_identical.pdf` |
| Supplement: identical fraction of sequence and gene content | `figures/figure_supp_identical_fraction_and_gene_content.py` | `bioapp_supp_identical_fraction.pdf` |
| Supplement: amino acid identity of orthologs of the examined pairs | `figures/figure_supp_pair_divergence.py` | `bioapp_supp_pair_divergence.pdf` |
| Supplement: *Anaerococcus* sp. Marseille-P3625 and *P. duerdenii* | `figures/figure_shared_regions.py pairs/peptoniphilaceae_2` | `bioapp_peptoniphilaceae_second_pair.pdf` |

Changing how the figures look:

- **Fonts, font sizes, figure width, colours:** `figures/plotstyle.py` (one place for all figures).
- **Panel sizes of the pair figures:** the constants at the top of `figures/figure_shared_regions.py`.
- **Which regions are drawn as gene maps, which genes are named, the names of the genomes:**
  `settings.ini` of the pair. **Labels and grouping of the compared genomes:** `genomes.tsv` of the pair.
- **Which pairs appear in the two supplementary figures of the examined pairs:** `pairs/examined_pairs.tsv`
  (delete a line to drop a pair; the pairs are labelled by name).

Each script writes a PDF for the manuscript and a PNG for a quick look.

## What is in the repository

```
analysis/        analysis scripts, numbered in the order in which they run; common.py holds shared helpers
figures/         one script per figure; plotstyle.py; make_all_figures.sh
pairs/           one folder per genome pair examined in detail
    <pair>/genomes.tsv      the two genomes of the pair (roles A and B) and the genomes compared with them
    <pair>/settings.ini     prefix of the pair's LaTeX macros, names, figure options
    <pair>/results/         result tables of the pair (our run)
    examined_pairs.tsv      the seven pairs, their identifier among the 140 pairs, and a number that names their macros in overview.tex
all_pairs/       the 140 genome pairs of the hierarchical edge bundling figure
    pairs.tsv               the pairs and their FracMinHash dN/dS in the manuscript's figure
    genomes.tsv             GTDB R214 metadata of their 243 genomes (taxonomy, isolate/MAG/SAG, CheckM estimates)
    make_inputs.py          re-creates pairs.tsv and genomes.tsv from their public sources
    results/                all_pairs_summary.tsv (one row per pair) and the tables of steps 13-15 (our run)
output/          everything the manuscript uses (our run)
    figures/                the figures (PDF and PNG)
    numbers/                LaTeX macros, one file per pair and overview.tex; every number quoted in the text
    data/                   supplementary data files S1 and S2
checks/          check_against_gtdb.py, check_expected_values.py, checksums of the GTDB R214 gene calls
environment.yml, reproduce_all.sh, export_to_manuscript.sh
```

Created by the scripts and not part of the repository (see `.gitignore`): `data/` (genomes, gene calls, reads),
`work/` (intermediate files) and `external/` (the FracMinHash dN/dS code).

## Analysis steps

For one pair examined in detail, `bash analysis/run_pair.sh pairs/<name>` runs steps 1 to 10.
The header of each script describes its method, inputs and outputs.

| Step | Script | What it does | Main output |
|---|---|---|---|
| 0 | `00_get_lab_code.sh` | clones `dnds-using-fmh` at the commit used in the manuscript | `external/` |
| 1 | `01_download_genomes.sh` | downloads the genomes of a pair from NCBI by assembly accession | `data/genomes/` |
| 2 | `02_call_genes.py` | Prodigal 2.6.3 gene calls with GTDB's options and translation tables | `data/genes/` |
| 3 | `03_shared_genes.py` | per gene, the fraction of its 21-mers shared with the other genome; near-identical genes | `genes_A.tsv`, `genes_B.tsv`, `genome_stats.tsv`, `contigs_A.tsv` |
| 4 | `04_fmh_dnds.py` | FracMinHash dN/dS with all genes, without near-identical genes (three cut-offs), and with every k-mer | `fmh_dnds.tsv` |
| 5 | `05_conventional_dnds.py` | alignment-based dN/dS (YN00) of reciprocal-best-hit orthologs | `orthologs.tsv` |
| 6 | `06_align_genomes.sh` | nucmer: genome A against genome B and the compared genomes | `work/<pair>/alignments/` |
| 7 | `07_shared_regions.py` | shared regions (E1, E2, ...), their genes and functional categories | `regions.tsv`, `region_genes.tsv`, `blocks.tsv`, `partner_contigs.tsv` |
| 8 | `08_scaffold_gaps.py` | scaffold gaps in or near the regions | `gaps_near_regions.tsv` |
| 9 | `09_presence_in_other_genomes.py` | presence of each region in genome B and the compared genomes | `presence.tsv` |
| 10 | `10_pair_numbers.py` | LaTeX macros of the pair | `output/numbers/<pair>.tex` |
| 11 | `11_all_pairs_run.sh` | steps 1-5 for each of the 140 pairs | `work/all_pairs/<pair>/results/` |
| 12 | `12_all_pairs_summary.py` | one row per pair; class of each pair | `all_pairs/results/all_pairs_summary.tsv` |
| 13 | `13_gene_removal_experiment.py` | FracMinHash dN/dS of conspecific genomes after random removal of genes | `all_pairs/results/gene_removal.tsv` |
| 14 | `14_translation_table_check.py` | the one pair whose genes GTDB called with translation table 4 | `all_pairs/results/translation_table_check.tsv` |
| 15 | `15_overview_numbers.py` | LaTeX macros about the 140 pairs; model of an identical fraction of sequence | `output/numbers/overview.tex`, `all_pairs/results/identical_fraction_model.tsv`, `carrier_genomes.tsv` |
| 16 | `16_map_source_reads.sh` | optional: maps a metagenome to the two MAGs of pair 1 | `work/.../reads/mapped.bam` |
| 17 | `17_read_level_summary.py` | optional: read depth and read pairs at the shared region | `read_depth.tsv`, `read_depth_windows.tsv`, `read_pairs.tsv` |
| 18 | `18_supplementary_data.py` | supplementary data files | `output/data/` |

### The optional read-level steps

Steps 16 and 17 ask which of the two organisms of pair 1 (*Methanobrevibacter* sp. RGIG2411 and
*Candidatus* Saccharibacteria bacterium RGIG2249) carries their shared region, by mapping the
metagenome from which RGIG2411 was assembled (SRA run SRR12529041, 106 million read pairs).
They need about 100 GB of disk and about one hour on 32 threads, and are therefore not part of
`reproduce_all.sh`. Their results are in the repository
(`pairs/methanobrevibacter_saccharibacteria/results/read_*.tsv`). To reproduce them:

```bash
export PYTHONPATH=analysis:figures
THREADS=32 bash analysis/16_map_source_reads.sh pairs/methanobrevibacter_saccharibacteria
python analysis/17_read_level_summary.py pairs/methanobrevibacter_saccharibacteria
python analysis/10_pair_numbers.py pairs/methanobrevibacter_saccharibacteria
python figures/figure_mag_misbinning.py
```

## Numbers in the manuscript

Every number that the manuscript quotes about these analyses is a LaTeX macro in
`output/numbers/*.tex`, written by steps 10 and 15, with a comment that says what it is:

```latex
\providecommand{\FRFmhDnds}{1.85}  % FracMinHash dN/dS of the pair, all genes
```

`bash export_to_manuscript.sh <manuscript folder>` copies the macros, the figures and the
supplementary data into the manuscript folder. The manuscript's section files list the value
each macro had when the text was written. Before exporting new numbers,

```bash
python checks/check_expected_values.py <manuscript folder>            # which quoted values changed?
python checks/check_expected_values.py <manuscript folder> --update   # accept the new values
```

shows which sentences need to be re-read.

## Data and versions

- **Genomes:** NCBI assemblies with version-pinned accessions (`pairs/*/genomes.tsv`,
  `all_pairs/genomes.tsv`; the complete list is `output/data/S2_Data_genomes.tsv`).
- **Consistency with the manuscript's data:** the manuscript's FracMinHash dN/dS values were
  computed from the gene calls of GTDB release 214. `checks/check_against_gtdb.py` confirms
  that the gene calls made here are identical to those of GTDB R214 for 242 of the 243 genomes
  (in the remaining genome the start codon of one gene differs by 9 bp) and that the FracMinHash
  dN/dS of all 140 pairs equals the value in the manuscript's figure (largest difference 2e-14).
- **FracMinHash dN/dS:** `KoslickiLab/dnds-using-fmh`, commit
  `ee94eb781bcd031eecfe87573fb7970d75fc5559`; k_aa = 7, k_nt = 21, scale factor 500, containment
  threshold 0.
- **Tools** (`environment.yml`): sourmash 4.9.4, sourmash_plugin_branchwater 0.9.15,
  Prodigal 2.6.3, MUMmer 4.0.1, DIAMOND 2.2.8, KaKs_Calculator 2.0.1, ncbi-datasets-cli 18.41.0,
  Python 3.12, Biopython 1.88, pandas 3.1.1, numpy 2.5.4, matplotlib 3.11.2; for the optional
  steps Bowtie2 2.5.5, SAMtools 1.25, sra-tools 3.4.1.
- **`all_pairs/pairs.tsv`** is the table behind the manuscript's hierarchical edge bundling
  figure (`figures/connect_species.csv` of
  [KoslickiLab/DnDs-visualization](https://github.com/KoslickiLab/DnDs-visualization), commit
  `060896e`), with a pair identifier added. **`all_pairs/genomes.tsv`** holds the GTDB R214 (214.1)
  metadata of its genomes. Both ship with the repository; `python all_pairs/make_inputs.py`
  re-creates them from these public sources (120 MB download).
- **`checks/gtdb_r214_gene_calls_md5.tsv`** holds checksums of the gene calls of the 243 genomes in
  the GTDB R214 files `gtdb_proteins_aa_reps_r214.tar.gz` and `gtdb_proteins_nt_reps_r214.tar.gz`
  (sequences joined, upper case, stop symbols removed from proteins).

## Definitions used throughout

- **Near-identical gene:** at least half of its nucleotide 21-mers occur in the coding sequences
  of the other genome, and it has at least 100 21-mers (about 97% nucleotide identity).
- **Fully shared contig:** at least 80% of its genes are near-identical.
- **Shared region:** nucmer alignments of at least 1 kb and 95% identity, merged when within
  10 kb of each other; named E1, E2, ... along genome A.
- **Present in another genome:** alignments of at least 90% identity cover at least 90% of the region.
- **Classes of genome pairs:** header of `analysis/12_all_pairs_summary.py`.

## Limitations

- Gene functions are NCBI product names sorted into broad categories by keyword
  (`analysis/07_shared_regions.py`); they colour the figures and are not a curated annotation.
- The direction of transfer is not determined for any pair.
- For *L. pacificus* and *S. pentaromativorans* only one assembly of the second genome exists and
  it has scaffold gaps at the shared regions.
