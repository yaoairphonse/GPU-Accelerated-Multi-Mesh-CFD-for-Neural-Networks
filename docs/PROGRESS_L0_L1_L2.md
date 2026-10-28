# Current progress map: L0, L1, and L2

Updated: 2026-10-28

This document is the short project-status map.  In this document, **L0/L1/L2
are readiness layers**, not mesh refinement indices:

- `LF1` and `LF2` remain the two fast 2-D airfoil mesh profiles.
- CRM AMR levels 9 and 10 remain Cartesian refinement levels.
- L0/L1/L2 below describe what the project can presently be trusted to do.

For a visual version, open
[the HTML quick guide](current_progress_quick_guide.html).

## One-page summary

| Layer | Purpose | Current status | Appropriate use today |
|---|---|---|---|
| **L0 — numerical trust** | Verify geometry, conservation, positivity, restart, and CPU/CUDA agreement | **Verified for the retained regressions** | Development gates, smoke tests, and detecting code regressions |
| **L1 — fast trend data** | Produce inexpensive 2-D or sectional CFD inputs for NNLCI and design screening | **Usable inside the tested envelope** | NACA4412-like low-Mach lift, moment, and surface-pressure trends with LF1/LF2 |
| **L2 — full 3-D aircraft** | Predict full-aircraft flow, surface pressure, and integrated loads on tetrahedral, hybrid, or Cartesian-AMR meshes | **Functional but not aerodynamically validated** | Mesh/solver research, GPU throughput, qualitative flow and Cp, and convergence development |

## L0 — numerical trust layer

L0 answers: **Does the discrete machinery behave consistently?**

Retained evidence:

- 3-D sphere conservation, positivity, restart, and CPU/CUDA regression.
- Uniform-flow preservation on the procedural CRM AMR mesh.
- Full-mesh geometry closure and global conservation audits.
- Fixed-step CPU/CUDA state agreement at approximately machine precision.
- Exact checkpoint reload and deterministic restart tests.
- Direct true-mesh visualization without scattered interpolation.

L0 is suitable as a software verification gate.  It does **not** establish
aircraft aerodynamic accuracy by itself.

## L1 — fast NNLCI trend-data layer

L1 answers: **Can the code quickly provide a useful aerodynamic trend?**

The strongest retained result is the 2-D Cartesian cut-cell NACA4412 case at
Mach 0.2:

| Profile | Controls | RTX 4090 process time | CL difference vs matched XFOIL | Cp RMS difference |
|---|---:|---:|---:|---:|
| LF1 | 4,718 | 20.38 s | -1.55% | 0.01282 |
| LF2 | 9,517 | 33.42 s | -1.77% | 0.01176 |

Both profiles passed the configured CL/CM/Cp stability gate.  They are useful
for rapid qualitative studies, NNLCI input generation, and version-to-version
comparisons within a similar tested regime.

The matrix-free block DP-LUR path runs on CPU and CUDA. At the conservative
CUDA CFL 0.5, LF1 reaches the gate in 5,600 pseudo-iterations and 3.22 s;
LF2 reaches it in 11,000 iterations and 6.25 s. Their CL/CD/CM values are
0.660704/0.006871/0.107404 and 0.660900/0.003068/0.106780. Uniform
preservation, perturbation damping, CPU/CUDA agreement, deterministic restart,
and face/load budget closure are regression tested. CFL 4 is not a portable
LF2 setting because it produced a non-decaying residual cycle.

Current limits:

- The comparison is against matched XFOIL, not wind-tunnel experiment.
- Drag and near-wall viscous accuracy are not validated.
- Broad Mach-number, geometry-family, and separated-flow validation remains.
- RAE2822 Cartesian and wedge/hexa paths prove mesh/solver capability, but the
  retained transonic results are not yet an accuracy claim.

## L2 — full 3-D aircraft layer

L2 answers: **Can the code ingest or generate a large 3-D aircraft mesh and
advance a compressible viscous solution on the GPU?**  The answer is yes, but
the current CRM loads are not converged validation data.

Implemented 3-D mesh paths:

- explicit-face tetrahedral, hybrid wedge/hexahedral, and general polyhedral;
- procedural multi-block Cartesian AMR with exact cut faces, agglomeration,
  2:1 interfaces, wall distance, restart, and CPU/CUDA backends;
- retained support for imported NASA tetrahedral CRM grids while Cartesian AMR
  remains the primary full-aircraft development path.

Current headline CRM AMR result:

| Item | Value |
|---|---:|
| Active controls | 21,523,036 |
| Minimum Cartesian spacing | 0.00390625 MAC |
| Latest continuation | 1,200 iterations; 11,900 cumulative |
| Latest continuation compute time | 549.22 s (9.15 min) |
| Saved endpoint CD/CL/CM | 0.029967 / 0.572087 / -0.041837 |
| Full low-to-low cycle mean CD/CL/CM | 0.028936 / 0.565420 / -0.052798 |
| NASA USM3D CD/CL/CM reference | 0.02995 / 0.54239 / -0.06249 |
| Endpoint global/near-body max residual | 3.76e-5 / 1.31e-3 |

The endpoint differs from NASA by CD/CL/CM =
+0.057%/+5.475%/+33.05%.  The full trajectory cycle differs by
-3.385%/+4.246%/+15.509%.  The solution therefore still fails the load-window
and fixed-point convergence requirements.  These are steady pseudo-iteration
diagnostics, not physical-time URANS averages.

Three recent opt-in experiments remain available but disabled by default:

1. The pressure-only normal probe is conservative and bounded, but changed
   wing Cp RMSE by only 0.018%; it is rejected as the primary accuracy fix.
2. The bounded SA face scheme reduced the SA residual but did not change the
   load trajectory materially.
3. The five-variable block DP-LUR neighbor action is verified, but added about
   10.4% compute cost for only 1--2% residual improvement in the matched test.

## What can be claimed now

- The repository is a working CPU/CUDA multi-mesh CFD research code.
- It handles 2-D Cartesian cut cells and multiple 3-D mesh representations.
- It demonstrates GPU operation from Mach 0.2 through Mach 1.65 on retained
  cases.
- L0 software/numerical regressions and the narrow L1 NACA4412 trend case are
  the strongest claims.
- L2 demonstrates full-aircraft meshing, throughput, flow visualization, and
  continued numerical development; it is not yet a validated CRM load solver.

## Immediate next work

Before another large CRM mesh or long run:

1. Use the implemented residual/load budgets to isolate procedural internal,
   same-level/2:1 interface, stored cut, wall, far-field, and SA-source effects.
2. Track upper-wing shock location and sectional pitching moment through the
   low-frequency load cycle.
3. Run matched lower-CFL or stronger-damping branches to distinguish slow
   pseudo-time damping from a numerical or genuinely unsteady mode.
4. Expand validation of the persistent 2-D CUDA DP-LUR path across broader
   airfoil families and transonic Mach numbers for high-throughput NNLCI input
   generation.
5. Only after the CRM diagnosis, reduce AMR memory and plan targeted 30--40M-cell
   refinement around the shock, leading/trailing edges, junction, tip, and
   wake.

## Where to look next

- [Current numerical status](STATUS.md)
- [Theory and operation manual](MANUAL.md)
- [Example matrix and true-mesh images](../examples/README.md)
- [CRM Cartesian-AMR runbook](../examples/3d_crm_amr/RUNBOOK.md)
- [Development history](DEVELOPMENT_HISTORY.md)
