# Native 2-D Cartesian-AMR solver

The native 2-D solver is implemented in the root [`src/`](../../src/) tree and
is built as `quick_airfoil`.  Its orchestration entry point is
[`quick_cfd.py`](../../quick_cfd.py).

Principal source files:

- `src/static_amr2d.cpp` and `.h`: mesh/control construction and CPU advance;
- `src/static_amr2d_cuda.cu`: CUDA residual and time-integration path;
- `src/static_amr2d_reconstruction.h`: limited primitive reconstruction;
- `src/static_amr2d_viscous.cpp`: viscous fluxes;
- `src/static_amr2d_sa.cpp`: SA-neg turbulence terms;
- `src/static_amr2d_boundary.cpp`: farfield/outlet boundary states;
- `src/static_amr2d_io.cpp`: fields, wall data, and aerodynamic loads.

Supported native 2-D features include LF1/LF2 quadtree refinement, exact
polygon cut cells and wall segments, conservative small-cell treatment,
CPU/CUDA execution, RK3 local pseudo-time, SA-neg or laminar flow, optional
low-Mach preconditioning, limited second-order mean-flow reconstruction, SI
exports, surface Cp, and integrated CL/CD/CM.

Start with the [NACA4412 runbook](../../examples/2d_naca4412/RUNBOOK.md) or the
[RAE2822 Cartesian-AMR example](../../examples/2d_rae2822_amr/).
