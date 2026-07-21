#!/bin/bash
set -euo pipefail

# Directory where config files will be written
config_dir="configs"
mkdir -p "${config_dir}"

# Constants
scaled=500
cds_input_list="dataset.csv"
translate_cds="no"
mode="bwpair"
cores=100

# Variable combinations
ks=(7 9 11)
thresholds=(0.1 0.2 0.02 0.05)

for k in "${ks[@]}"; do
    for threshold in "${thresholds[@]}"; do

        sample="gtdb_protein_rep_pairwise_${scaled}scale_${threshold}threshold_k${k}"
        wd="/data/${sample}"
        out="${sample}_out"

        config_file="${config_dir}/${sample}.conf"

        cat > "${config_file}" <<EOF
# Automatically generated configuration file

k=${k}
scaled=${scaled}
threshold=${threshold}

sample=${sample}
wd=${wd}
cds_input_list=${cds_input_list}

translate_cds=${translate_cds}
mode=${mode}
cores=${cores}

out=${out}
EOF

        echo "Created ${config_file}"

    done
done
