Here, you can reproduce FMH Omega estimations. Produce you config files by modifying the handy dandy `config_generation.sh` file to fit your environment and parameters.

Generate config files
```bash config_generation.sh```

Before reproducing estimations, make sure that you have the following:

A `dataset.csv` file:

This file includes the file pathways of genome and protein fasta files for each genome.

FMH Omega Program available on your environment:

Please refer to our FMH Omega repo here: [FMH Omega](https://github.com/KoslickiLab/dnds-using-fmh/tree/dev_jzr2/src)

Run config files as desired using the following command
```./run_fmh_omega.sh ${config file name}```

