# Development status — 2026-10-24

## Retained results

| Case | Size | Result | Status |
|---|---:|---|---|
| NACA4412 LF1 | 4,718 CV | CL 0.66444, Cp RMS 0.01282 vs XFOIL | gate passed |
| NACA4412 LF2 | 9,517 CV | CL 0.66297, Cp RMS 0.01176 vs XFOIL | gate passed |
| 3-D sphere | 7,680 cells | CPU/CUDA about 1e-15 | regression passed |
| RAE2822 hybrid LF1 | 5,051 cells | CL/CD/CM 0.803917/0.096155/-0.098244 | solver gate passed; reference consolidation pending |
| CRM tetra | 10,067,380 cells | CL/CD/CM 0.581806/0.052443/-0.143849 | unconverged |
| CRM AMR base | 7,661,606 controls | CL/CD/CM 0.335620/0.025291/+0.131069 | unconverged |
| CRM AMR wall-refined | 14,011,654 controls | CL/CD/CM 0.357993/0.028121/-0.049050 | unconverged |
| CRM AMR upper-wing L10 | 21,523,036 controls | late-1,000 pseudo-iteration average CL/CD/CM 0.523661/0.026297/-0.024508 | 5,000 cumulative updates; unconverged |

NASA CRM USM3D comparison values are CL/CD/CM
0.54239/0.02995/-0.06249.  The displayed CRM AMR result is a trapezoidal
average over cumulative iterations 4,000--5,000, not a selected endpoint and
not a physical-time URANS average.  The last-500 average is
0.514545/0.025687/-0.012119, so load-window convergence has not been reached.
The saved endpoint is 0.511187/0.026104/-0.017449 and remains a raw transient
sample.  Wing Cp RMSE at that saved state is 0.100340.

## Current blockers

1. Low-frequency CRM load/shock motion, especially pitching moment; distinguish
   slow pseudo-time damping from a numerical or physical unsteady mode.
2. AMR memory per control before a 30--40M mesh.
3. Full second-order viscous/SA verification.
4. Long-window load and pressure/viscous drag convergence.
5. External experimental/Fluent reference packaging.

## Next development phase

Three isolated options are now implemented and disabled by default: a
pressure-only wall-normal probe, a bounded positive/dissipative SA face
operator, and a matrix-free five-variable Euler neighbor action for procedural
DP-LUR.  CPU/CUDA and conservation tests pass.  Matched CRM tests rejected the
pressure probe and block action as stand-alone accuracy fixes; the bounded SA
operator lowers SA residuals but does not remove the load mode.  Full decision
records are under `docs/development_history/`.

Before another arbitrary long run, instrument residual and load budgets by
procedural internal, 2:1 interface, stored cut, wall, and far-field face
families.  Track upper-wing shock position and sectional moment, then use
matched lower-CFL/stronger-damping branches to determine whether the observed
motion is a pseudo-time artifact.  A physical-time URANS average requires a
fixed physical time step and converged inner iterations; steady pseudo-
iterations must not be presented as physical-time samples.  Memory reduction
and a targeted 30--40M refinement remain later work.
