# Simple airfoil interface

`simple_airfoil.py` is the minimal entry point for design sweeps and NNLCI
input generation.  The user supplies the geometry plus only six CFD choices:

1. algorithm (`rk3` or `dplur`);
2. nominal mesh size (4,700 or 9,500 controls);
3. Mach number;
4. angle of attack in degrees;
5. freestream static pressure in Pa;
6. freestream static temperature in K.

The output directory is also required.  The script automatically detects an
`x,y` GXXX CSV or an ordered `x,y,patch` polygon CSV.  It locks spatial order,
SA-neg, Sutherland viscosity, Turkel low-Mach preconditioning, CFL, relaxation
sweeps, convergence limits, and backend to retained presets.

## Command

```bash
python3 simple_airfoil.py \
  --geometry examples/2d_naca4412/geometry.csv \
  --algorithm rk3 \
  --mesh-cells 4700 \
  --mach 0.2 \
  --aoa 1.5 \
  --pressure-pa 101325 \
  --temperature-k 288.15 \
  --out output/simple_naca4412
```

`rk3` selects the retained CUDA path.  `dplur` selects the CPU reference
path.  The nominal 4,700/9,500 selections map to LF1/LF2; the actual cut-cell
count depends slightly on geometry and is always reported.  Geometry is
normalized by a fixed one-metre reference chord in this simple interface.
The underlying solver also supports CUDA DP-LUR through
`config/naca4412_m02_order2_dplur_cuda.json` or the offline case builder; it is
kept outside this deliberately minimal two-choice wrapper so the DP-LUR option
continues to work on machines without CUDA.

The main machine-readable result is `simple_summary.json`.  Normal field,
surface, load, mesh, convergence-history, and HTML outputs remain available in
the same output directory.

## Measured runtime

Fresh runs on the project workstation and RTX 4090 used NACA4412, Mach 0.2,
AOA 1.5 degrees, 101325 Pa, and 288.15 K.  All four combinations passed the
same raw-residual plus stable CL/CM/Cp gate.

| Algorithm | Nominal / actual controls | Iterations | Process time | CL | CD | CM (CCW-positive) |
|---|---:|---:|---:|---:|---:|---:|
| CPU DP-LUR | 4,700 / 4,718 | 1,900 | 7.70 s | 0.660046 | 0.006801 | 0.107250 |
| CUDA RK3 | 4,700 / 4,718 | 43,427 | 20.38 s | 0.664440 | 0.005959 | 0.107868 |
| CUDA RK3 | 9,500 / 9,517 | 70,800 | 33.42 s | 0.662969 | 0.003748 | 0.107167 |
| CPU DP-LUR | 9,500 / 9,517 | 5,900 | 48.25 s | 0.658080 | 0.003898 | 0.106460 |

For the coarse preset, CPU DP-LUR is currently fastest.  For the fine preset,
CUDA RK3 is currently fastest because LF2 requires a lower locked DP-LUR CFL
to avoid the residual oscillation observed with the LF1 setting.  Do not use
iteration counts to compare algorithms directly: one DP-LUR outer update
contains five implicit relaxation sweeps, while RK3 has three explicit stages.

The timing table is hardware- and geometry-specific; every production run
records its own process and pipeline times.  Lift, moment, and surface pressure
are the supported trend quantities.  Drag is exported but is not yet an
aerodynamic validation claim for this rapid cut-cell workflow.
