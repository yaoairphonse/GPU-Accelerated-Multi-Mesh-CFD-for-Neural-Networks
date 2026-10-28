#!/usr/bin/env bash
set -euo pipefail

executable="${1:-build_cpu/quick_airfoil}"
python3 quick_cfd.py run \
  --points examples/2d_naca4412/geometry.csv --point-format polygon \
  --config examples/2d_naca4412/case_dplur.json \
  --out output/2d_naca4412_dplur --levels LF1 \
  --backend cpu --executable "${executable}"

echo "Open output/2d_naca4412_dplur/index.html"
