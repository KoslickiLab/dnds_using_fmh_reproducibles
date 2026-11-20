#!/usr/bin/env bash
# Draw every figure of the biological applications from the result tables in this repository.
# No genomes need to be downloaded for this. Output: output/figures/*.pdf (and .png previews).
#
#   usage:  bash figures/make_all_figures.sh
#
# To change the look of the figures (fonts, sizes, colours, width), edit figures/plotstyle.py.
set -euo pipefail
export PYTHONPATH=analysis:figures

python figures/figure_mag_misbinning.py                                   # bioapp_mag
python figures/figure_shared_regions.py pairs/faecalibacterium_roseburia  # bioapp_faecalibacterium_roseburia
python figures/figure_shared_regions.py pairs/bacteroidales               # bioapp_bacteroidales
python figures/figure_shared_regions.py pairs/peptoniphilaceae_1          # bioapp_peptoniphilaceae
python figures/figure_shared_regions.py pairs/peptoniphilaceae_2          # bioapp_peptoniphilaceae_second_pair (supplement)
python figures/figure_shared_regions.py pairs/roseburia                   # bioapp_roseburia
python figures/figure_shared_regions.py pairs/rhodobacteraceae            # bioapp_rhodobacteraceae
python figures/figure_supp_identical_fraction_and_gene_content.py         # bioapp_supp_identical_fraction (supplement)
python figures/figure_supp_pair_divergence.py                             # bioapp_supp_pair_divergence (supplement)
python figures/figure_supp_near_identical_genes.py                        # bioapp_supp_near_identical (supplement)
