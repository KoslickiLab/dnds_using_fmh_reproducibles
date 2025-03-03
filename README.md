# This repo is for the reproducibles of manuscript on fmh dnds

## Simulation dataset generation:

Generate random sequences:


Apply mutation on real sequences:


From the simulated random sequences, produce an AXT from FASTA file and apply KaKs_Calculator on generated AXT dataset.

    python helper_scripts/mutation_simulation_on_random_data/apply_kaks_calculator_on_random_sequence.py

## Figures:

Figure 2A

    python helper_scripts/mutation_simulation_on_random_data/pairplot_against_NG86_different_lengths_and_ksizes_v3.py

Figure 2B

    python helper_scripts/mutation_simulations_on_real_data/pairplot_against_NG86_different_ksizes_LAMA3.py

Figure 3A

## Supplemental Figures:

Figure 6:

    python helper_scripts/mutation_simulation_on_random_data/histogram_kaks_test.py

Figure 7:

    python helper_scripts/mutation_simulations_on_real_data/histogram_different_scales_ecoli.py