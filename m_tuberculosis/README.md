# FracMinHash dN/dS across a collection of *Mycobacterium tuberculosis* isolates

Code, result tables, figures and LaTeX numbers for the supplementary vignette on closely related
genomes of the manuscript *Leveraging FracMinHash Containment for Genomic dN/dS*.

Gan et al. (2025, mSystems 10:e00898-25) report an ordering of the groups of genes of
*M. tuberculosis* (essential genes below non-essential genes below the PE/PPE family) and restrict
their pairwise dN/dS analysis to 5,000 random isolates because of its cost. This folder computes
FracMinHash dN/dS for every pair of 39,379 of their isolates (775,333,131 pairs; 49,904 isolates and
1,245,179,656 pairs including the read-derived ones) for the whole coding genome and for each group of
genes, compares it with `seqinr`'s `kaks` (the estimator Gan et al. used) on the same sequences, and
measures how sketching (scale factors 2 to 1,000) degrades the values when isolates differ at about
one site in several thousand. Only broad evolutionary patterns are examined; no gene or isolate is
characterised individually.

Every FracMinHash dN/dS value comes from
[KoslickiLab/dnds-using-fmh](https://github.com/KoslickiLab/dnds-using-fmh) at commit `19a25de`
(pull request 18), which adds `--mode reference`: for many genomes that are all close to one
reference, it computes the containments of every pair exactly from each genome's differences to the
reference sketch, instead of comparing every pair of sketches in full. The manuscript's estimator,
parameters and output columns are unchanged.

## Reproduce everything

```bash
conda env create -f environment.yml
conda activate fmh_dnds_mtb
bash reproduce_all.sh
```

This is a computation for a server: it downloads 16 GB of AllTheBacteria assembly archives (165 GB
once extracted), needs about 5 TB of disk for `work/` and ran on 250 cores (`CORES=250`, step 32).
The tool's own run times are in `results/trends_tool/timings.tsv` and `results/tool/sketching.tsv`.
All scripts are run from the folder that contains this README.

The results of our run are part of the repository (`results/`), so the figures and the numbers can be
regenerated without the computation:

```bash
python figures/42_figures.py     # the eight panels, from the tables in results/
python analysis/43_numbers.py    # results/numbers/mtb.tex
```

## Figures

Each panel is its own file (`results/figures/panels/`, PDF and PNG), drawn by `figures/42_figures.py`.
A panel in which a legend or a label covers plotted data is not written (`figures/plotstyle.py`,
`save_panel`). The panels are put together, with their letters, by the LaTeX in
`latex/mtb_figures.tex`; captions are in `latex/mtb_captions.tex`.

| Figure | Panels | What they show |
|---|---|---|
| FracMinHash dN/dS over all pairs, every k-mer | `mtb_main_a` to `mtb_main_d` | whole coding genome by relatedness of the pair; the groups of genes beside Gan et al.'s values; against `seqinr` on the same sequences; against the divergence of the pair |
| What sketching costs | `mtb_scale_a` to `mtb_scale_d` | containments of draft assemblies against anchored sequences; rank correlation with every k-mer, medians and the orderings of Gan et al. at scale factors 1 to 1,000 |

## What is in the repository

```
analysis/        analysis scripts, numbered in the order in which they run (gaps are steps of the
                 investigation that the vignette does not use); common.sh holds shared settings
figures/         42_figures.py (every panel) and plotstyle.py
checks/          31a-31b: reference mode against main mode; 46: every panel recomputed from the raw output
results/         result tables of our run (see "Result tables"), the panels and the number file
latex/           the vignette's text, captions, figure layout and bibliography (standalone.tex
                 compiles them on their own)
environment.yml, environment_reads.yml, reproduce_all.sh, export_to_manuscript.sh
```

Created by the scripts and not part of the repository (see `.gitignore`): `data/`, `work/`,
`external/` (the FracMinHash dN/dS code), `logs/` and `tmp/`.

## Analysis steps

The header of each script describes its method, inputs and outputs.

| Step | Script | What it does | Main output |
|---|---|---|---|
| 0 | `analysis/00_get_tools.sh` | clones `dnds-using-fmh` at the commit with `--mode reference` | `external/` |
| 1 | `analysis/01_get_atb_metadata.sh` | AllTheBacteria release 0.2 metadata, checked against its MD5 sums | `data/atb/` |
| 1b | `analysis/01b_get_published_data.sh` | supplements of Gan et al. 2025 and DeJesus et al. 2017, the MTBC ancestral sequence, H37Rv (GCF_000195955.2); MD5-checked | `data/gan2025/`, `data/essential/`, `data/reference/` |
| 2 | `analysis/02_match_gan_to_atb.py` | the runs of Gan et al. (Table S1) matched to AllTheBacteria assemblies | `data/isolates/gan_to_atb.tsv` |
| 3 | `analysis/03_get_atb_assemblies.sh` | the assembly archives that hold those isolates | `data/atb/0.2/assembly/` |
| 5 | `analysis/05_ena_runinfo.py` | ENA metadata of every run | `data/isolates/ena_runinfo.tsv` |
| 6 | `analysis/06_read_runs.py` | optional: runs for the read pipeline (runs without an assembly, and 300 calibration runs) | `data/isolates/read_pipeline_runs.txt` |
| 7, 8 | `analysis/07_reads_to_snvs.sh`, `08_run_read_pipeline.sh` | optional: substitutions from reads (bwa, bcftools; thresholds of Gan et al.) | `work/snv_reads/` |
| 10 | `analysis/10_extract_assemblies.sh` | extracts the assemblies of the isolates | `data/assemblies/` |
| 11 | `analysis/11_call_genes.sh` | Prodigal gene calls of the draft assemblies (used by step 12 only) | `work/genes/` |
| 12 | `analysis/12_smoke_scales.py` | 50 draft assemblies compared through their own gene calls, with the tool's unchanged main mode | `results/smoke/` |
| 13 | `analysis/13_gene_groups.py` | essential, non-essential and PE/PPE genes and the PE/PPE high-quality segments, from published tables only | `results/groups/` |
| 14, 14b | `analysis/14_variants_from_assemblies.sh`, `14b_coverage.py` | substitutions of every assembly against the ancestral sequence (minimap2, paftools) | `work/variants/`, `results/inputs/coverage.tsv` |
| 24 | `analysis/24_site_ratio.py` | ratio of synonymous to non-synonymous sites of the reference coding sequences | `results/site_ratio.tsv` |
| 30 | `analysis/30_write_inputs.py` | coding sequences and proteins of every isolate, per group of genes | `work/fasta/`, `results/inputs/` |
| 31a, 31, 31b | `checks/31a_check_inputs.py`, `31_check_reference_mode.sh`, `31b_check_summary.py` | main mode and reference mode on the same isolates (1,000 and 80 chosen at random): every pair, every column identical | `results/checks/` |
| 32 | `analysis/32_run_reference_mode.sh` | FracMinHash dN/dS of every pair, whole coding genome and each group, scale factors 1 to 1,000 | `work/tool/` |
| 33 | `analysis/33_conventional.py` | `seqinr`'s `kaks` on the same sequences, 200 random isolates (19,900 pairs) | `results/conventional_tool/` |
| 34 | `analysis/34_reads_vs_assembly.py` | optional: isolates with both reads and an assembly, the two inputs compared | `results/reads_vs_assembly_tool/` |
| 35 | `analysis/35_sketch_times.py` | sketching times, from the log of step 32 | `results/tool/sketching.tsv` |
| 40 | `analysis/40_summarise_tool_output.py` | distributions per run and scale factor; a subsample of pairs shared by all runs | `results/tool/` |
| 40b | `analysis/40b_exact_quantiles.py` | exact quantiles, overall, per category of pair and per bin of divergence | `results/tool/*/quantiles*.tsv` |
| 41 | `analysis/41_trends.py` | the comparisons of the vignette | `results/trends_tool/` |
| 44 | `analysis/44_divergence_tail.py` | the most divergent pairs and the isolates behind them | `results/trends_tool/divergence_tail*.tsv` |
| 42 | `figures/42_figures.py` | the panels | `results/figures/panels/` |
| 43 | `analysis/43_numbers.py` | every number the text quotes, as LaTeX macros | `results/numbers/mtb.tex` |
| 46 | `checks/46_review_figures.py` | every plotted quantity recomputed from the raw output | `results/review/checks.tsv` |

### The optional read pipeline

10,525 of the isolates have no assembly in AllTheBacteria; their substitutions were called from the
reads (steps 06-08: SRA download, bwa mem to the ancestral sequence, bcftools; mapping quality > 30,
depth > 10, allele frequency > 0.9, as in Gan et al.). They enter only the run over every isolate
(`all_cds_every_isolate`), which shows that including them barely moves the median whole-genome
value; every other analysis uses the assembly-derived isolates. The pipeline downloads the reads
of about 10,900 runs and is not part of `reproduce_all.sh`. Without it, step 30 writes the
assembly-derived isolates only and the run over every isolate covers them only; our results for it are
in `results/tool/all_cds_every_isolate/` and `results/reads_vs_assembly_tool/`. To run it:

```bash
conda env create -f environment_reads.yml
conda activate fmh_dnds_mtb_reads
python analysis/06_read_runs.py          # (in the main environment)
JOBS=80 bash analysis/08_run_read_pipeline.sh
conda activate fmh_dnds_mtb
python analysis/34_reads_vs_assembly.py  # after step 32
```

## Numbers in the manuscript

Every number the vignette quotes is a LaTeX macro in `results/numbers/mtb.tex`, written by step 43
from the result tables, with a comment that says what it is:

```latex
% median dN/dS, whole coding genome, every pair of the assembly-derived isolates
\providecommand{\mtbMedianWhole}{...}
```

`bash export_to_manuscript.sh <manuscript folder>` copies the macros, the panels and the vignette's
LaTeX into the manuscript folder. To build the vignette on its own:

```bash
mkdir -p latex/numbers latex/figures
cp results/numbers/mtb.tex latex/numbers/ && cp results/figures/panels/mtb_*.pdf latex/figures/
cd latex && tectonic standalone.tex      # or pdflatex + bibtex
```

## Checks

- **Reference mode equals main mode** (step 31, `results/checks/reference_mode.tsv`): on 1,000
  isolates (PE/PPE high-quality segments; scale factors 1, 2, 5) and 80 isolates (whole coding genome;
  scale factors 1, 5), plus the reference, both modes report identical values for every pair and every
  column (1,507,980 pairs over the five comparisons). Main mode cannot be run on more whole coding
  genomes at scale factor 1: its single sketch file passes 4 GB, and sourmash cannot read the ZIP64
  archive it then becomes.
- **Every panel against the raw output** (step 46, `results/review/checks.tsv`): quantiles found by a
  different algorithm, pair-by-pair comparisons on a different subsample, `seqinr` re-run on
  two-sequence alignments, and every value drawn from a table looked up in the tool's per-pair output.

## Data and versions

- **Isolates:** the 49,954 runs (49,949 distinct) of Gan et al. 2025, Table S1; assemblies from
  AllTheBacteria release 0.2.
- **Reference:** the inferred MTBC ancestral sequence (Comas et al. 2013) as distributed with the code
  of Gan et al. ([irongump/MTB_PEPPE](https://github.com/irongump/MTB_PEPPE), commit `4cf8101`), which
  has the coordinates of H37Rv; coding sequences from the H37Rv RefSeq annotation (GCF_000195955.2).
- **Groups of genes:** essentiality from DeJesus et al. 2017 (Table S3); PE/PPE membership from the
  H37Rv annotation; high-quality PE/PPE segments from Gan et al. 2025 (Tables S2, S3).
- **FracMinHash dN/dS:** `KoslickiLab/dnds-using-fmh`, commit
  `19a25de8d214e6797983e3dfdcdf958e7ec9ddc1`; k_aa = 7 (k_nt = 21), max containment, constant
  0.77/2.23, threshold 0, scale factor 1 unless stated.
- **Tools** (`environment.yml`): sourmash 4.9.4, sourmash_plugin_branchwater 0.9.15, minimap2 2.31,
  Prodigal 2.6.3, R 4.5.3 with seqinr 4.2-44, Python 3.12, numpy 2.5.5, pandas 3.1.1, scipy 1.18.1,
  pyarrow 26.0.1, matplotlib 3.11.2; for the read pipeline (`environment_reads.yml`) sra-tools 3.4.1,
  seqtk 1.5, bwa 0.7.19, samtools 1.25, bcftools 1.25.

## Limitations

- Each isolate is the ancestral sequence with its substitutions applied: insertions, deletions and
  genes absent from the reference are not represented, and positions an assembly does not cover keep
  the reference base.
- The absolute values depend on the estimator and its constant; they are compared with each other, and
  with `seqinr` on the same sequences, not with the per-gene pN/pS values of Gan et al., which come
  from a different estimator.
