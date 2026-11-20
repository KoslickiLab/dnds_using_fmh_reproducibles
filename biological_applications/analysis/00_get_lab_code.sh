#!/usr/bin/env bash
# Step 0. Fetch the FracMinHash dN/dS implementation at the exact commit used in the manuscript.
# All analyses call it unmodified (external/dnds-using-fmh/fmh_dnds.py).
#
#   usage:  bash analysis/00_get_lab_code.sh
set -euo pipefail
commit=ee94eb781bcd031eecfe87573fb7970d75fc5559
if [ ! -d external/dnds-using-fmh ]; then
    git clone https://github.com/KoslickiLab/dnds-using-fmh.git external/dnds-using-fmh
fi
git -C external/dnds-using-fmh checkout --quiet $commit
echo "FracMinHash dN/dS code at commit $(git -C external/dnds-using-fmh rev-parse --short HEAD)"
