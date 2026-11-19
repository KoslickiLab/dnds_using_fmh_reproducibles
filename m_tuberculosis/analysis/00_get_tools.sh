#!/usr/bin/env bash
# Step 0. The FracMinHash dN/dS implementation, at the commit that adds --mode reference
# (KoslickiLab/dnds-using-fmh, pull request 18). Every FracMinHash dN/dS value of this vignette comes
# from external/dnds-using-fmh/fmh_dnds.py.
#
#   usage:  bash analysis/00_get_tools.sh
set -euo pipefail
commit=19a25de8d214e6797983e3dfdcdf958e7ec9ddc1
if [ ! -d external/dnds-using-fmh ]; then
    git clone https://github.com/KoslickiLab/dnds-using-fmh.git external/dnds-using-fmh
fi
# the commit is on the default branch once pull request 18 is merged, and on its branch until then
git -C external/dnds-using-fmh fetch --quiet origin "$commit" 2>/dev/null \
    || git -C external/dnds-using-fmh fetch --quiet origin pull/18/head
git -C external/dnds-using-fmh checkout --quiet "$commit"
echo "FracMinHash dN/dS code at commit $(git -C external/dnds-using-fmh rev-parse --short HEAD)"
