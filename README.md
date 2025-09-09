# This repo is for the reproducibles of manuscript on fmh dnds

# Table of Contents

- [Environment setupe](#Environment-Setup)
    - [Conda](#Conda-environment)
    - [Pip](#Pip-environment)
- [Datasets](#Datasets)
    - [Simulations](#Simulations)
    - [GTDB](#GTDB)
- [Reproduce dN/dS estimations](#Reproduce)
    - [FracMinHash dN/dS](#fracminhash-dnds-estimations)
    - [Traditional dN/dS models](#traditional-dnds-models)

# Environment setup

Create either a conda or pip environment to reproduce this figure.

## Conda environment

```
conda env create -f environment.yml
conda activate reproduce
```

## Pip environment

```
python3 -m venv reproduce_pip_env
source reproduce_pip_env/bin/activate
pip install -r requirements.txt
```

# Datasets 

## Simulations

While using a scale factor of 1, we produce FracMinHash dN/dS estimation employing varying k-sizes to compare with the traditional dn/ds model, NG86. Please follow instructions here to generate these and more simulations.

## GTDB 

<!-- ### Random

To generate random sequence simulations please refer to:


#### Negative

#### Positive

### Real (LAMA3)

#### Negative

#### Positive

## Real dataset (GTDB)
-->

# Reproduce dN/dS estimations

## FracMinHash dN/dS



## Traditional dN/dS models

From the simulated sequences, produce an AXT from FASTA file and apply KaKs_Calculator [1] on generated AXT dataset.

```
python helper_scripts/mutation_simulation_on_random_data/apply_kaks_calculator_on_random_sequence.py
```

# Main figures

These are the scripts to generate the main figures of our manuscript.

## Figure 2

Figure 2A represents how well FracMinHash dN/dS estimations are being made when compared to the traditional dN/dS model, NG86, on random sequences. Additionally, the figure compares varying k-sizes and sequence lengths. Please execute the following command to produce the figure for random sequence simulations.

```
python figure_scripts/figure2a.py
```

![Figure 2A](manuscript_figures/figure2a.png)

Similar to Figure 2A, Figure 3b represents how well FracMinHash dN/dS estimations are being made when compared to the traditional dN/dS model, NG86, but this time we run simulations on a real sequence. The figure also compares varying k-sizes. Please execute the following command to produce the figure for a real sequence simulations.

```
python figure_scripts/figure2b.py
```

![Figure 2B](manuscript_figures/figure2b.png)

## Figure 3

Figure 3A

    python figure_scripts/disk_usage_figure.py

Figure 3B

    python figure_scripts/runtimes_stackplot_figure.py

## Figure 4

## Figure 5

To reproduce Figure 5A, the hierarchical edge bundling figure, please refer to the following repo: [Hierarchichal Edge Bundling](https://github.com/KoslickiLab/DnDs-visualization)

![Figure 5A: Hierarchical Edge Bundling Figure](https://github.com/KoslickiLab/DnDs-visualization/blob/main/figures/output_species.png?raw=true)


# Supplemental figures

Figure 6:

    python helper_scripts/mutation_simulation_on_random_data/histogram_kaks_test.py

Figure 7:

    python helper_scripts/mutation_simulations_on_real_data/histogram_different_scales_ecoli.py

# References

[1] Zhang, Z., Li, J., Zhao, X.-Q., Wang, J., Wong, G. K.-S., & Yu, J. (2006). KaKs_Calculator: Calculating Ka and Ks through model selection and model averaging. Genomics, Proteomics & Bioinformatics, 4(4), 259–263. Oxford University Press.

# Please cite

Leveraging FracMinHash Containment for Genomic dN /dS. Judith S. Rodriguez, Mahmudur Rahman Hera, and David Koslicki. In preparation.