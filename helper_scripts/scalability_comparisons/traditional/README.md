# Testing scalibility across dnds models

## Traditional method

### Sample 
python sampling.py

We start to record time from the alignments

### alignments

```
nohup /usr/bin/time -v python3 alignment.py   --input /data/jzr5814/dnds_scalability_comparison/sample_5/fasta   --out_dir /data/jzr5814/dnds_scalability_comparison/iter1/sample_5/alignments > /data/jzr5814/dnds_scalability_comparison/sample_5/alignment.log 2>&1 & 
```

```
nohup /usr/bin/time -v python3 alignment.py   --input /data/jzr5814/dnds_scalability_comparison/sample_10/fasta   --out_dir /data/jzr5814/dnds_scalability_comparison/iter1/sample_10/alignments > /data/jzr5814/dnds_scalability_comparison/sample_10/alignment.log 2>&1 &
```

```
nohup /usr/bin/time -v python3 alignment.py   --input /data/jzr5814/dnds_scalability_comparison/sample_100/fasta   --out_dir /data/jzr5814/dnds_scalability_comparison/iter1/sample_100/alignments > /data/jzr5814/dnds_scalability_comparison/sample_100/alignment.log 2>&1 &
```

```
nohup /usr/bin/time -v python3 alignment.py   --input /data/jzr5814/dnds_scalability_comparison/sample_1000/fasta   --out_dir /data/jzr5814/dnds_scalability_comparison/iter1/sample_1000/alignments > /data/jzr5814/dnds_scalability_comparison/sample_1000/alignment.log 2>&1 &
```

### covert alignments to axt files and concatenate axt files

```
nohup /usr/bin/time -v python3 axt_convert.py \
  --input  /data/jzr5814/dnds_scalability_comparison/sample_5/alignments \
  --out_dir /data/jzr5814/dnds_scalability_comparison/sample_5/axt \
  >> /data/jzr5814/dnds_scalability_comparison/sample_5/axt.log 2>&1 &
```

```
nohup /usr/bin/time -v python3 axt_convert.py \
  --input  /data/jzr5814/dnds_scalability_comparison/sample_10/alignments \
  --out_dir /data/jzr5814/dnds_scalability_comparison/sample_10/axt \
  >> /data/jzr5814/dnds_scalability_comparison/sample_10/axt.log 2>&1 &
```

```
nohup /usr/bin/time -v python3 axt_convert.py \
  --input  /data/jzr5814/dnds_scalability_comparison/sample_100/alignments \
  --out_dir /data/jzr5814/dnds_scalability_comparison/sample_100/axt \
  >> /data/jzr5814/dnds_scalability_comparison/sample_100/axt.log 2>&1 &
```

```
nohup /usr/bin/time -v python3 axt_convert.py \
  --input  /data/jzr5814/dnds_scalability_comparison/sample_1000/alignments \
  --out_dir /data/jzr5814/dnds_scalability_comparison/sample_1000/axt \
  >> /data/jzr5814/dnds_scalability_comparison/sample_1000/axt.log 2>&1 &
```

### run kaks calculator

```
nohup /usr/bin/time -v python3 kaks_calculator_job.py   --base_dir /data/jzr5814/dnds_scalability_comparison   --samples 5,10,100,1000   --methods NG,YN   --input_name "axt/sequences.axt"   --out_template "sequences_{method}.axt.kaks"   --log_template "sequences_{method}.kaks.log"   >> /data/jzr5814/dnds_scalability_comparison/kaks_grid.log 2>&1 &
```

### Summary statistics 

Please find here: traditiona_dnds_walltime_summary.xlsx