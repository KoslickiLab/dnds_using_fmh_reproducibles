# Generate simulated datasets for FracMinHash dN/dS validation

# Table of Contents

- [Environment setup](#environment-setup)
- [Datasets](#datasets)

# Environment setup

Before simulating selected sequences, you'll need to set up an environment. Please find instructions here: [Environment setup](https://github.com/KoslickiLab/dnds_using_fmh_reproducibles/tree/main#Environment-Setup)

# Datasets 

In our manuscript, we used a scale factor of 1 and varying k-sizes to validate that FracMinHash dN/dS is accurately estimating selection and to compare with the traditional dn/ds model results. Please find the datasets we have used for our manuscript here: XYZ

# Simulate

While using a scale factor of 1, we produce FracMinHash dN/dS estimation employing varying k-sizes to compare with the traditional dn/ds model, NG86. Please follow instructions here to generate these and more simulations.

```
python random_selection_simulation.py --len 5000 --prate 0.01 --wd ../
python random_selection_simulation.py --len 10002 --prate 0.01 --wd ../
python random_selection_simulation.py --len 20001 --prate 0.01 --wd ../
```