# NACA4412 LF1/LF2 runbook

## Purpose

This is the retained fast 2-D accuracy example.  It exercises polygon input,
Cartesian cut cells, small-cell treatment, low-Mach preconditioning, SA-neg,
limited second-order mean-flow reconstruction, LF1/LF2 execution, surface Cp,
and the CL/CM/Cp convergence gate.

Conditions: Mach 0.2, alpha 1.5 degrees, Re 4.6595 million, one-metre chord.

## 1. Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCFD_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=89
cmake --build build -j 8 --target quick_airfoil
```

Use `-DCFD_CUDA=OFF -DCFD_OPENMP=ON` and `build_cpu` for a CPU build.

## 2. Inspect inputs

- `geometry.csv` is the ordered polygon with patch information.
- `case.json` requests second-order reconstruction, Turkel preconditioning,
  fixed CFL behavior, SA-neg, and the CL/CM/Cp convergence metric.
- `reference/expected_result.json` records the retained XFOIL comparison.

## 3. Run both fidelity levels

```bash
bash examples/2d_naca4412/run.sh cuda build/quick_airfoil
```

The Python driver validates and resolves the configuration, constructs LF1 and
LF2 meshes, launches the standalone solver, rejects stale binaries, and writes
an HTML result index.  Expected control-volume counts are about 4,718 and
9,517.  Fresh simple-interface runs required 20.38 and 33.42 process seconds.

## 4. Check the result

```bash
python3 validation/compare_xfoil_surface.py \
  --cp external_data/xfoil_naca4412_ref/cp.txt \
  --polar external_data/xfoil_naca4412_ref/polar.txt \
  --run output/2d_naca4412 \
  --out output/2d_naca4412/xfoil_comparison.json \
  --plot output/2d_naca4412/xfoil_cp_comparison.png
```

The exact comparison CLI depends on the retained XFOIL reference layout; run
`python3 validation/compare_xfoil_surface.py --help` before substituting a new
reference.  The recorded LF1/LF2 lift errors are -1.55%/-1.77%, moment-
magnitude errors +1.19%/+0.53%, and Cp RMS errors 0.01282/0.01176.

## 5. Convergence interpretation

The configured gate requires the raw residual, CL/CM change, and wall-length-
weighted Cp change to remain below their thresholds for five reports.  Passing
this development gate does not validate skin-friction drag or transition.

## 6. CPU DP-LUR reference

The opt-in 2-D DP-LUR implementation is currently a CPU reference path.  It
uses the same mesh, residual, second-order reconstruction, Turkel
preconditioning, SA-neg model, and convergence gate:

```bash
bash examples/2d_naca4412/run_dplur.sh build_cpu/quick_airfoil
```

The retained LF1 run reached the gate after 1,900 pseudo-iterations in 7.70 s
process time, with CL/CD/CM = 0.660046/0.006801/0.107250 and no rejected
updates.  LF2 uses a lower locked CFL and reached the gate after 5,900
iterations in 48.25 s.  Lift and moment are trend checks; drag remains
unvalidated.  Setting
`output.audit_budgets=true` also writes per-face-family conservation and
chordwise pressure/viscous load budgets next to the internal flow files.

CUDA DP-LUR is not implemented for this 2-D path yet.  Use `case.json` and
`run.sh` for the retained CUDA RK3 production path.

## 7. CPU/CUDA comparison

Run the same short fixed-step case with both backends and compare the exported
states.  The retained 19-step order-2 test agreed to relative Linf 1.20e-15.
