# Development Status — 2026-10-24

## Benchmark & Verification Summary

| Case | Grid size | Recorded result | Status |
|---|---:|---|---|
| CRM AMR upper-wing L10 | 21,523,036 active cells | 4,000–5,000 iter mean $C_L/C_D/C_M = 0.523661 / 0.026297 / -0.024508$, wing $C_p$ RMSE 0.10034 | 5,000 cumulative updates; low-frequency $C_M$ mode active |
| CRM AMR wall-refined | 14,011,654 active cells | $C_L/C_D/C_M = 0.357993 / 0.028121 / -0.049050$ | Short-window checkpoint |
| CRM AMR base | 7,661,606 active cells | $C_L/C_D/C_M = 0.335620 / 0.025291 / +0.131069$ | Coarse baseline checkpoint |
| CRM tetrahedral | 10,067,380 tetrahedra | $C_L/C_D/C_M = 0.581806 / 0.052443 / -0.143849$ | Unconverged throughput run |
| 3D sphere ($M = 1.65$) | 7,680 polyhedra | CPU/CUDA agreement within $\sim 10^{-15}$ | Regression passed |
| RAE 2822 hybrid LF1 | 5,051 wedge/hex cells | $C_L/C_D/C_M = 0.803917 / 0.096155 / -0.098244$ | Solver verification passed |
| NACA 4412 LF1 | 4,718 cut cells | $C_L = 0.66444$, $C_p$ RMS 0.01282 vs. XFOIL | Validated |
| NACA 4412 LF2 | 9,517 cut cells | $C_L = 0.66297$, $C_p$ RMS 0.01176 vs. XFOIL | Validated |

For the 21.52M-cell NASA CRM run at Mach 0.85, $\alpha = 2.5^\circ$, $\text{Re}_{\text{MAC}} = 5\times 10^6$, the NASA USM3D reference coefficients are $C_L = 0.54239$, $C_D = 0.02995$, and $C_M = -0.06249$. Trapezoidal averaging over cumulative pseudo-time iterations 4,000–5,000 yields $C_L/C_D/C_M = 0.523661 / 0.026297 / -0.024508$, whereas the 4,500–5,000 window gives $0.514545 / 0.025687 / -0.012119$ (final iteration 5,000 state: $0.511187 / 0.026104 / -0.017449$, wing surface $C_p$ RMSE $0.100340$). Because these are local pseudo-time iterations and the late window still shows pitching-moment drift, load convergence has not yet been reached.

## Current Technical Focus

1. Damping the low-frequency CRM pitching-moment and shock-position oscillation in pseudo-time.
2. Reducing per-cell GPU memory footprint prior to 30–40M cell AMR grids.
3. Completing full second-order viscous and SA turbulence verification.
4. Long-window pressure and viscous drag stabilization across multiple mesh refinement levels.

## Recent Solver Evaluations

Three optional solver extensions were implemented and verified for conservation and CPU/CUDA parity (all disabled by default):
- **Wall-normal pressure probe (`--wall-state-reconstruction normal-probe`):** Precomputes two fluid-side donors for bounded one-sided pressure extrapolation. On a matched 21.52M-cell CRM A/B test, it altered wing $C_p$ RMSE by only 0.018% and was left off by default.
- **Bounded SA face scheme (`--sa-face-scheme bounded`):** Combines neighbor-bounded convection, positive bounded non-orthogonal diffusion, and a consistent Rusanov/two-point scalar DP-LUR update, reducing the 300-step SA residual by $1.77\times$.
- **5×5 block Euler neighbor coupling (`--dplur-block-coupling 1`):** Applies a matrix-free 5-variable Euler normal-flux Jacobian action during procedural DP-LUR sweeps (CPU/CUDA parity $1.87\times 10^{-16}$). On the 21.52M-cell CRM grid it reduced residuals by 1–2% at a 10.4% compute overhead.

Detailed evaluation logs are stored under `docs/development_history/`.
