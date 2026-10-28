# AMR solver methods and verification update — 2026-10-28

This note records the numerical changes added after the retained 5,000-step
CRM run and the corresponding 2-D implementation intended for fast NNLCI
input generation.  It distinguishes verified code behavior from aerodynamic
validation.

## Current method map

| Component | 2-D Cartesian cut cell | 3-D Cartesian octree AMR |
|---|---|---|
| Spatial residual | Shared subface finite-volume flux, exact cut-wall segments, order 1 or limited order 2 | Procedural block interiors, stored cut faces, exact wall fragments, and explicit 2:1 interfaces |
| Steady marching | GPU/CPU SSP-RK3; new CPU matrix-free block DP-LUR reference | CUDA/CPU DP-LUR, five sweeps and scalar or five-variable Euler neighbor action |
| Low-Mach treatment | Optional Turkel residual preconditioning | Not used for the Mach 0.85 CRM case |
| Viscous/turbulence model | Laminar or SA-neg with Sutherland viscosity | SA-neg with wall-distance table and optional bounded SA face operator |
| Wall pressure | Cut-segment pressure; optional bounded two-donor normal probe | Cut-fragment pressure; optional pressure-only normal probe |
| Diagnostics | Face-family residual budgets and 24-bin chordwise pressure/viscous load budgets | Block, same-level, 2:1, cut, far-field, wall, and SA-source budgets plus 24-bin x/y load budgets |

The normal-probe wall-pressure option is deliberately disabled by default.
It is a bounded diagnostic reconstruction, not a wall function, and matched
CRM tests did not identify it as an accuracy fix.

## New 2-D DP-LUR reference path

The 2-D path uses local pseudo-time, repeated matrix-free relaxation sweeps,
a scalar spectral/viscous backbone, and a four-variable Euler
normal-flux Jacobian action when block coupling is enabled.  SA remains a
scalar implicit backbone.  A rejected non-positive update halves the working
CFL.  The saved state contains the pseudo-iteration counter, so a 10+10-step
restart is deterministic and matches a continuous 20-step run exactly in the
unit test.

The current verified implementation is CPU only.  The existing CUDA 2-D
production path remains SSP-RK3; a CUDA DP-LUR port has not yet been
implemented or benchmarked.

Run the retained LF1 example from the repository root:

```bash
cmake -S . -B build_cpu -DCMAKE_BUILD_TYPE=Release \
  -DCFD_CUDA=OFF -DCFD_OPENMP=ON
cmake --build build_cpu -j 8
bash examples/2d_naca4412/run_dplur.sh build_cpu/quick_airfoil
```

The configuration lives in `examples/2d_naca4412/case_dplur.json`.  The key
controls are:

```json
{
  "solver": {
    "integrator": "dplur",
    "dplur": {
      "CFL": 4.0,
      "sweeps": 5,
      "omega": 0.8,
      "block_coupling": true,
      "maximum_iterations": 5000,
      "check_interval": 100
    }
  },
  "output": {"audit_budgets": true}
}
```

On the retained NACA 4412 LF1 case, the CPU reference reached the existing
CL/CM/Cp stability gate after 1,900 pseudo-iterations in 7.70 s process time,
with no rejected update.  A separate fresh LF2 run required 5,900 iterations
and 48.25 s with its lower locked CFL:

| Controls | Iterations | Time | CL | CD | CM (CCW-positive) |
|---:|---:|---:|---:|---:|---:|
| 4,718 | 1,900 | 7.70 s | 0.660046 | 0.006801 | 0.107250 |
| 9,517 | 5,900 | 48.25 s | 0.658080 | 0.003898 | 0.106460 |

The matched XFOIL values are CL 0.6749, CD 0.00572, and moment magnitude
0.1066.  Lift and moment are suitable narrow-envelope trend checks.  The
cut-cell viscous drag is not yet validated, so the CD agreement must not be
treated as an accuracy claim.

Unit verification covers exact uniform-flow preservation, damping of a
pressure perturbation, and exact 20-step versus 10+10-step restart equality.
For a short audited airfoil state, same-level and 2:1 signed internal-face
integrals cancel to machine zero, and the binned pressure/viscous loads sum to
the integrated result within approximately 1e-14.

## CRM continuation and face-family diagnostics

The 21,523,036-control CRM mesh was continued from 5,000 to 10,700 updates at
CFL 4, followed by 1,200 updates at CFL 12.  The last segment required 549.22
s on one RTX 4090, or about 0.458 s per outer update.  Its saved endpoint is:

| Total pseudo-iterations | CD | CL | CM | Global max residual | Near-body max residual |
|---:|---:|---:|---:|---:|---:|
| 11,900 | 0.029967 | 0.572087 | -0.041837 | 3.76e-5 | 1.31e-3 |

NASA USM3D reference values are CD/CL/CM =
0.02995/0.54239/-0.06249.  Endpoint differences are +0.057%, +5.475%, and
33.05% respectively when normalized by the signed reference coefficient.
The endpoint is not proof of convergence: the full observed low-to-low load
cycle, cumulative iterations 6,550--10,200, averages
CD/CL/CM = 0.028936/0.565420/-0.052798.  This is a steady
pseudo-iteration diagnostic average, not a physical-time URANS average.

The 3-D executable now accepts:

- `--face-budget-audit 1` to write `face_budget.csv` and
  `wall_load_budget.csv`;
- `--write-state 0|1` to suppress or retain the large checkpoint;
- `--write-surface 0|1` to suppress or retain triangle Cp/load output.

The face audit executes one family at a time, allowing same-level and 2:1
interfaces to be checked independently.  The audit-on and audit-off solution
states were bitwise identical in the verification run.

## Interpretation and remaining work

The new diagnostics show that the late CRM trajectory is much closer to the
reference than the 5,000-step state, especially for drag and lift.  Pitching
moment and the near-body residual remain insufficiently stable.  The current
result therefore supports the procedural-AMR method as a throughput and trend
development path, but it is not a validated steady CRM solution.

Before increasing the mesh to 30--40 million controls, use the face and
sectional load budgets to localize the low-frequency mode.  If that mode
cannot be damped as a steady iteration, add a true physical-time dual-time
path and report a physical-time average.  CUDA DP-LUR for the 2-D NNLCI path
is a separate performance task; it should preserve the CPU reference results
before becoming the default.
