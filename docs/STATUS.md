# Development Status — 2026-10-28

For the short readiness map, see the [L0/L1/L2 progress summary](PROGRESS_L0_L1_L2.md) or open the [visual HTML quick guide](current_progress_quick_guide.html).

## Benchmark & Verification Summary

| Case | Grid size | Recorded result | Status |
|---|---:|---|---|
| CRM AMR upper-wing L10 | 21,523,036 active cells | endpoint $C_D/C_L/C_M = 0.029967 / 0.572087 / -0.041837$ | 11,900 cumulative updates; low-frequency $C_M$ mode active |
| CRM AMR wall-refined | 14,011,654 active cells | $C_L/C_D/C_M = 0.357993 / 0.028121 / -0.049050$ | Short-window checkpoint |
| CRM AMR base | 7,661,606 active cells | $C_L/C_D/C_M = 0.335620 / 0.025291 / +0.131069$ | Coarse baseline checkpoint |
| CRM tetrahedral | 10,067,380 tetrahedra | $C_L/C_D/C_M = 0.581806 / 0.052443 / -0.143849$ | Unconverged throughput run |
| 3D sphere ($M = 1.65$) | 7,680 polyhedra | CPU/CUDA agreement within $\sim 10^{-15}$ | Regression passed |
| RAE 2822 hybrid LF1 | 5,051 wedge/hex cells | $C_L/C_D/C_M = 0.803917 / 0.096155 / -0.098244$ | Solver verification passed |
| NACA 4412 LF1 | 4,718 cut cells | $C_L = 0.66444$, $C_p$ RMS 0.01282 vs. XFOIL | Validated |
| NACA 4412 LF2 | 9,517 cut cells | $C_L = 0.66297$, $C_p$ RMS 0.01176 vs. XFOIL | Validated |
| NACA 4412 CUDA DP-LUR LF1/LF2 | 4,718 / 9,517 cut cells | 3.22 / 6.25 s process time | Both residual + CL/CM/Cp gates passed; drag not validated |
| NACA 4412 CPU DP-LUR LF1/LF2 | 4,718 / 9,517 cut cells | 7.70 / 48.25 s process time | Portable reference; both trend gates passed |

For the 21.52M-cell NASA CRM run at Mach 0.85, $\alpha = 2.5^\circ$, $\text{Re}_{\text{MAC}} = 5\times 10^6$, the NASA USM3D reference is $C_D/C_L/C_M = 0.02995/0.54239/-0.06249$. The 11,900-step endpoint differs by +0.057%/+5.475%/+33.05%. A full observed low-to-low trajectory cycle over iterations 6,550–10,200 averages $0.028936/0.565420/-0.052798$, or -3.385%/+4.246%/+15.509%. This is a steady pseudo-iteration diagnostic, not a physical-time URANS average; global/near-body maximum residuals remain $3.76\times10^{-5}/1.31\times10^{-3}$.

## Current Technical Focus

1. Damping the low-frequency CRM pitching-moment and shock-position oscillation in pseudo-time.
2. Reducing per-cell GPU memory footprint prior to 30–40M cell AMR grids.
3. Completing full second-order viscous and SA turbulence verification.
4. Long-window pressure and viscous drag stabilization across multiple mesh refinement levels.
5. Broader Mach/geometry validation of the new CUDA DP-LUR 2D NNLCI path.

## Recent Solver Evaluations

Three optional solver extensions were implemented and verified for conservation and CPU/CUDA parity (all disabled by default):
- **Wall-normal pressure probe (`--wall-state-reconstruction normal-probe`):** Precomputes two fluid-side donors for bounded one-sided pressure extrapolation. On a matched 21.52M-cell CRM A/B test, it altered wing $C_p$ RMSE by only 0.018% and was left off by default.
- **Bounded SA face scheme (`--sa-face-scheme bounded`):** Combines neighbor-bounded convection, positive bounded non-orthogonal diffusion, and a consistent Rusanov/two-point scalar DP-LUR update, reducing the 300-step SA residual by $1.77\times$.
- **5×5 block Euler neighbor coupling (`--dplur-block-coupling 1`):** Applies a matrix-free 5-variable Euler normal-flux Jacobian action during procedural DP-LUR sweeps (CPU/CUDA parity $1.87\times 10^{-16}$). On the 21.52M-cell CRM grid it reduced residuals by 1–2% at a 10.4% compute overhead.

Detailed evaluation logs are stored under `docs/development_history/`.

Residual and load budgets are now available by procedural block, same-level/2:1 interface, stored cut, wall, far-field, and SA-source families. The 2D CPU/CUDA DP-LUR path adds matching deterministic-restart and face/load audit coverage. CUDA DP-LUR passes the retained NACA4412 LF1/LF2 gate at CFL 0.5; broad validation remains pending.
