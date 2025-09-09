# This repo is for the reproducibles of manuscript on fmh dnds

# Table of Contents

# Environment setup

# Simulation dataset generation:

Generate random sequences:


Apply mutation on real sequences:


From the simulated random sequences, produce an AXT from FASTA file and apply KaKs_Calculator on generated AXT dataset.

    python helper_scripts/mutation_simulation_on_random_data/apply_kaks_calculator_on_random_sequence.py

# Main figures

## Figure 2

## Figure 2 (Panel A)

    python helper_scripts/mutation_simulation_on_random_data/pairplot_against_NG86_different_lengths_and_ksizes_v3.py

![Figure 2A](dnds_using_fmh_reproducibles/manuscript_figures/figure2a.png)

Figure 2B

    python helper_scripts/mutation_simulations_on_real_data/pairplot_against_NG86_different_ksizes_LAMA3.py

Figure 3A

    python figure_scripts/disk_usage_figure.py

Figure 3B

    python figure_scripts/runtimes_stackplot_figure.py

## Figure 5

### Figure 5 (Panel A): Hierarchical Edge Bundling Figure

To reproduce hierarchical edge bundling figure, please refer to the following repo: [Hierarchichal Edge Bundling](https://github.com/KoslickiLab/DnDs-visualization)

![Figure 5A: Hierarchical Edge Bundling Figure](https://github.com/KoslickiLab/DnDs-visualization/blob/main/figures/output_species.png?raw=true)


# Supplemental figures

Figure 6:

    python helper_scripts/mutation_simulation_on_random_data/histogram_kaks_test.py

Figure 7:

    python helper_scripts/mutation_simulations_on_real_data/histogram_different_scales_ecoli.py
