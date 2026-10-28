# Quick CFD: GPU-Accelerated Multi-Fidelity Compressible Flow Solver

**2D/3D compressible RANS/Euler solver on Cartesian AMR cut-cell, tetrahedral, hybrid, and polyhedral meshes, with CPU/OpenMP and CUDA backends (RK3 and DP-LUR) for aerodynamic load and surface pressure prediction.**

Quick CFD is a GPU-focused solver built to generate multi-fidelity compressible flow fields and aerodynamic trends for neural-network fidelity-lift and surrogate models (such as NNLCI). It includes a 3D procedural octree-AMR cut-cell solver for full aircraft configurations, an unstructured 3D explicit-face solver for general polyhedral and hybrid grids, and a 2D Cartesian cut-cell solver for airfoils. The repository provides the C++/CUDA and Python sources, mesh converters, example runbooks, and a theory manual.

**Verification and validation status:** Geometry closure, freestream preservation, density/pressure positivity, checkpoint restarts, and fixed-step CPU/CUDA bitwise or near-machine-precision parity are verified across the test suite. The tested 2D NACA 4412 case is benchmarked against XFOIL for lift, moment, and surface-pressure trends. On the 3D NASA Common Research Model (CRM) 21.52M-cell octree-AMR grid, the complete observed low-to-low steady pseudo-iteration cycle (iterations 6,550--10,200) averages $C_D/C_L/C_M = 0.028936/0.565420/-0.052798$ against NASA USM3D $0.02995/0.54239/-0.06249$, corresponding to differences of $-3.385\%/+4.246\%/+15.509\%$. The instantaneous $C_M$ difference was $10.960\%$ at iteration 10,900, while the later unconverged 11,900-step endpoint drifted to $C_D/C_L/C_M = 0.029967/0.572087/-0.041837$. These are steady pseudo-iteration diagnostics rather than a physical-time average; pitching moment and the near-body residual remain unstable, so fixed-point convergence and quantitative 3D load validation are not claimed.

![CRM full-aircraft 3D surface pressure from the Cartesian-AMR solver](examples/3d_crm_amr/showcase/crm_amr_surface_cp_3d.png)

## 3D Full-Aircraft Showcase (NASA CRM, 21.52M Cells, Mach 0.85)

All aircraft surface and volume slices below are plotted directly on exact cut-wall triangles and procedural AMR control volumes without `scipy.interpolate.griddata` smoothing.

| Multi-view surface $C_p$ & load evolution | Horizontal wing slice ($z/\text{MAC} = 0.80$, $C_p$ & Mach) |
|---|---|
| ![CRM surface Cp projections and load history](examples/3d_crm_amr/showcase/crm_surface_cp.png) | ![CRM horizontal x-y wing Mach and Cp](examples/3d_crm_amr/showcase/crm_xy_wing_mach_cp.png) |
| **Full-span crossflow slice ($x/\text{MAC} = 5.0$, $C_p$ & Mach)** | **Near-symmetry slice ($y \to 0^+$, $C_p$ & Mach)** |
| ![CRM full-span y-z Mach and Cp](examples/3d_crm_amr/showcase/crm_yz_fullwing_mach_cp.png) | ![CRM symmetry x-z Mach and Cp](examples/3d_crm_amr/showcase/crm_xz_symmetry_mach_cp.png) |
| **3D octree-AMR block & cut-cell cross-sections** | **Spanwise AMR leaf blocks & wing cut-cell zoom** |
| ![CRM true AMR block and cell cross-sections](examples/3d_crm_amr/showcase/crm_amr_true_mesh_sections.png) | ![CRM true leaf blocks and individual AMR cells](examples/3d_crm_amr/showcase/crm_amr_yz_x5_mesh.png) |
| **Cumulative load & residual history** | **Exact-triangle 3D surface $C_p$** |
| ![CRM 11,900-step pseudo-iteration history](examples/3d_crm_amr/showcase/crm_amr_convergence_11900.png) | ![CRM three-dimensional Cp perspectives](examples/3d_crm_amr/showcase/crm_amr_surface_cp_3d.png) |

## Benchmarks and Throughput

GPU timings below were measured on one NVIDIA RTX 4090; the explicitly marked 2D CPU reference uses the same workstation. Values are end-to-end process or solver wall-clock times over test cases from Mach 0.2 to 1.65:

| Case | Mach | Mesh size | Recorded runtime | Status |
|---|---:|---:|---:|---|
| NASA CRM Cartesian AMR | 0.85 | 21,523,036 active cells (from 31.48M octree cells) | latest 1,200 CFL-12 updates in 549.22 s; 11,900 cumulative | Full-cycle mean errors: $-3.39\%/+4.25\%/+15.51\%$ in $C_D/C_L/C_M$; fixed point remains unconverged |
| NASA CRM tetrahedral | 0.85 | 10,067,380 tetrahedra | 5,850 DP-LUR sweeps in 600.20 s (610.35 s total) | Unstructured 3D throughput test |
| 3D sphere regression | 1.65 | 7,680 polyhedral cells | 2.15 s total process time | Conservation and CPU/CUDA parity check ($10^{-15}$) |
| RAE 2822 hybrid LF1 | 0.78 | 5,051 wedge/hex cells | 49,000 DP-LUR sweeps in 11.24 s | Solver and CPU/CUDA parity verified |
| NACA 4412 LF1 / LF2, CUDA RK3 | 0.20 | 4,718 / 9,517 cells | 20.38 / 33.42 s total process time | $C_L$, $C_M$, and surface $C_p$ match XFOIL reference |
| NACA 4412 LF1 / LF2, CPU DP-LUR reference | 0.20 | 4,718 / 9,517 cells | 7.70 / 48.25 s total process time | Both trend gates passed; 2D CUDA DP-LUR is not yet implemented |

On the 21.52M-cell CRM AMR grid, the latest five-sweep block-DP-LUR segment averaged about 0.458 s per outer update. The full observed low-to-low load cycle averages $C_D/C_L/C_M = 0.028936/0.565420/-0.052798$. This is a steady pseudo-iteration diagnostic rather than a physical-time URANS average. Building the full octree cut-cell mesh from the closed aircraft surface takes 8 min 52 s from scratch, or 5 min 20 s when refining from a cached parent grid.

## Mesh Representations

| Dim | Grid type | Application in this repo |
|---|---|---|
| 3D | Procedural multi-block octree AMR cut cells | $8\times 8\times 8$ leaf blocks with exact surface clipping, cell agglomeration, and 2:1 refinement interfaces |
| 3D | Explicit owner/neighbor polyhedral faces | Unstructured tetrahedra, prism/hex hybrid slabs, and 6-block cubed-sphere grids |
| 2D | Structured Cartesian & quadtree AMR cut cells | Multi-fidelity (LF1/LF2) airfoils with exact boundary segments and small-cell merging |

For the 3D octree AMR solver, regular interior block connectivity is evaluated procedurally on the fly; only block interfaces, physical boundaries, cut-face exceptions, and surface triangles are stored explicitly in device memory.

## Solver Paths

All solvers use a density-based finite-volume formulation (explicit RK3 or implicit DP-LUR) with low-Mach preconditioning for subsonic cases; pressure-based segregated schemes (SIMPLE/PISO) are not used.

| Path | Mesh topology | Discretization & time marching | Verification / validation |
|---|---|---|---|
| 3D CRM AMR | Multi-block octree cut cell | Procedural DP-LUR (scalar or $5\times 5$ block Euler action) | 7.66M, 14.01M, and 21.52M control-volume runs; 11,900-step state remains unconverged |
| 3D CRM tetra | Body-fitted unstructured tetrahedra | Explicit-face DP-LUR | 10.07M-cell import, checkpoint restart, gradient & SA tests |
| 3D sphere | 6-block cubed-sphere polyhedra | Explicit-face CPU/CUDA | Conservation, positivity, CPU/CUDA parity |
| 3D RAE 2822 | Mixed wedge/hex thin wing section | Explicit-face DP-LUR | LF1 run and CPU/CUDA parity |
| 2D airfoil | Cartesian cut cell (LF1/LF2) | CUDA/CPU RK3; CPU matrix-free block DP-LUR; SA-neg and low-Mach preconditioning | RK3 validated against XFOIL ($C_L, C_M, C_p$); DP-LUR regressions and LF1/LF2 trend gates passed |

## Build and Quick Start

For the locked no-tuning airfoil interface, choose only the algorithm, nominal mesh count, and freestream conditions:

```bash
python3 simple_airfoil.py \
  --geometry examples/2d_naca4412/geometry.csv \
  --algorithm rk3 --mesh-cells 4700 \
  --mach 0.2 --aoa 1.5 \
  --pressure-pa 101325 --temperature-k 288.15 \
  --out output/simple_naca4412
```

See the [simple interface and measured runtime table](docs/SIMPLE_AIRFOIL_INTERFACE.md).

CUDA build (defaults to `sm_89` for RTX 4090):

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCFD_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=89
cmake --build build -j 8
ctest --test-dir build --output-on-failure
```

CPU-only build (OpenMP):

```bash
cmake -S . -B build_cpu -DCMAKE_BUILD_TYPE=Release \
  -DCFD_CUDA=OFF -DCFD_OPENMP=ON
cmake --build build_cpu -j 8
ctest --test-dir build_cpu --output-on-failure
```

Run the self-contained 3D supersonic sphere test case:

```bash
bash examples/3d_sphere/run_smoke.sh build_cpu cpu
```

## Additional 3D Unstructured & 2D Airfoil Gallery

Body-fitted 3D unstructured grids and 2D airfoil test cases included in the verification suite (full index in the [example showcase](examples/README.md#visual-showcase)):

| Body-fitted 3D CRM tetrahedral section | 3D cubed-sphere polyhedra ($M = 1.65$) |
|---|---|
| ![CRM exact tetrahedron-plane intersections](examples/3d_crm_tetra/showcase/crm_tetra_exact_section.png) | ![Six-block sphere surface, exact plane cells, and mesh](examples/3d_sphere/showcase/sphere_m165_exact_mesh.png) |
| **Wedge/hex hybrid wing slab (RAE 2822)** | **2D quadtree AMR cut cells (NACA 4412)** |
| ![RAE2822 wedge and hexahedron hybrid grid](examples/3d_rae2822_hybrid/showcase/rae2822_lf2_hybrid.png) | ![NACA4412 true AMR members, flow, and wall Cp](examples/2d_naca4412/showcase/naca4412_lf2.png) |

## Documentation

- [L0/L1/L2 current-progress quick guide](docs/current_progress_quick_guide.html) and [text summary](docs/PROGRESS_L0_L1_L2.md)
- [2026-10-28 AMR methods, 2D DP-LUR, and CRM continuation update](docs/AMR_METHODS_2026-10-28.md)
- [No-tuning airfoil interface and runtime table](docs/SIMPLE_AIRFOIL_INTERFACE.md)
- [Theory and operation manual](docs/MANUAL.md)
- [Development log and history](docs/DEVELOPMENT_HISTORY.md)
- [Example matrix and runbooks](examples/README.md)
- [Mesh and flow visual showcase](examples/README.md#visual-showcase)
- [Current verification & convergence status](docs/STATUS.md)
- [Binary and mesh file formats](docs/FORMATS.md)
- [References and external datasets](docs/REFERENCES.md)

## Repository Layout

- `solvers/` — [2D/3D solver executables and source map](solvers/README.md)
- `solvers/3d/` — Procedural-AMR `amr3d` and explicit-face `sphere3d` solvers, mesh importers, restart tools, and unit tests
- `tools/cartesian_amr/` — Octree builder, exact triangle-box clipping, cut-cell moments, cell agglomeration, and GPU table packing
- `tools/hybrid/` — 2D mixed tri/quad extrusion to 3D wedge/hex hybrid grids
- `src/`, `pipeline/`, `quick_cfd.py` — full 2D `quick_airfoil` Cartesian cut-cell workflow
- `simple_airfoil.py` — locked six-input airfoil/NNLCI interface
- `examples/` — Reproducible runbooks and configuration files for each test case
- `validation/` — Post-processing and reference comparison scripts

## External Meshes and License

Raw NASA CRM and RAE 2822 archive grids are not stored in the repository; the scripts under `examples/` document the download URLs, checksums, and conversion commands. See [LICENSE_NOT_SELECTED.md](LICENSE_NOT_SELECTED.md) regarding current licensing status.
