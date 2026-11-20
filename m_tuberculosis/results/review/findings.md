# Adversarial review of the figures (2026-11-18/19)

Two parts: reading each panel and caption against the data and the source paper, and recomputing
every plotted quantity from the tool's raw per-pair output with independent code
(`scripts/46_review_figures.py`; results in `checks.tsv`, counts in `verdicts.json`).

## What the review changed

1. **Quantiles were bin centres.** Boxes, medians, quartiles and divergence bins came from histograms
   with bins 0.005 wide, so a three-decimal median could be off by up to 0.0025 (whole coding genome:
   0.567 drawn and quoted, 0.570 exact). Step 40b now computes exact quantiles; panels and numbers use
   them. A consequence: with the read-derived isolates the median moves from 0.570 to 0.568, so the text
   no longer says it is "unchanged".
2. **Gan et al.'s values were mislabelled.** 0.69, 0.84 and 0.90 are per-gene pN/pS counted from
   mutation events on phylogenies (Gan et al. 2025, Results), not pairwise `seqinr` values, and the paper
   does not call them medians. The text, the caption of main b and its legend now say so, and that only
   their ordering is comparable.
3. **Panel 1a legend.** "within a lineage" could be read as including pairs within a sublineage; the
   labels now name the three categories exactly, and the caption gives the bin width.
4. **Captions.** The anchored pairs of scale a are drawn from the pairs of scale b (not "from all
   pairs"); main c states that pairs identical over a group have no value with either estimator;
   numbers that were typed by hand (100,000 pairs per divergence bin, 4,000 points, 200,000 pairs and
   p < 10^-300 for the Wilcoxon tests) are now macros from script output.
5. **The review's own first run had two false failures**: it counted pairs with an undefined value
   (0/0, identical sequences) as missing from the raw output. Fixed; those pairs are undefined in both
   the table and the raw output, and `seqinr` gives no value for them either.

## What was recomputed (all pass)

- main a: every pair counted once; share below 1; the histogram drawn holds every finite value.
- main b: the five quantiles of every box, every group's pair count, the medians quoted in the text, and
  the three diamonds against Gan et al.'s values.
- main c: every FracMinHash value in the comparison tables against the raw output; seqinr dN/dS = ka/ks;
  `seqinr` re-run on 60 random pairs as two-sequence alignments; the pairs beyond the axes.
- main d: three divergence bins (lowest, middle, highest drawn) recomputed exactly.
- scale a: the 4,000 anchored points against the raw output; the draft-assembly points against the
  tool's own output of step 12.
- scale b and d: rank correlations and orderings recomputed on a different subsample of 12.1 million
  pairs (largest difference from the drawn value 0.0006, i.e. sampling error).
- scale c: the median of every group at every scale factor, exactly.

Not drawn: 0.019% of pairs fall outside the x axis of main a (note in `checks.tsv`).
