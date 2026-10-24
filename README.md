# Quick CFD: GPU-Accelerated Multi-Fidelity Compressible Flow Solver

**2D/3D compressible RANS/Euler solver on Cartesian AMR cut-cell, tetrahedral, hybrid, and polyhedral meshes, with CPU/OpenMP and CUDA backends (RK3 and DP-LUR) for aerodynamic load and surface pressure prediction.**

Quick CFD is a GPU-focused solver built to generate multi-fidelity compressible flow fields and aerodynamic trends for neural-network fidelity-lift and surrogate models (such as NNLCI). It includes a 2D Cartesian cut-cell solver for airfoils, an unstructured 3D explicit-face solver for general polyhedral and hybrid grids, and a 3D procedural octree-AMR cut-cell solver for full aircraft configurations. The repository provides the C++/CUDA and Python sources, mesh converters, example runbooks, and a theory manual.

**Verification and validation status:** Geometry closure, freestream preservation, density/pressure positivity, checkpoint restarts, and fixed-step CPU/CUDA bitwise or near-machine-precision parity are verified across the test suite. The 2D NACA 4412 case is validated against XFOIL. In 3D, the NASA Common Research Model (CRM) cases currently serve as large-grid GPU throughput benchmarks; long-window iterative convergence, SA turbulence closure tuning, and integrated drag validation ($C_L, C_D, C_M$) on the 3D CRM grids are still in progress.

![CRM full-aircraft 3D surface pressure from the Cartesian-AMR solver](examples/3d_crm_amr/showcase/crm_amr_surface_cp_3d.png)

## Benchmarks and GPU Throughput (RTX 4090)

Timings below were measured on a single NVIDIA RTX 4090 and reflect end-to-end process or solver wall-clock time over test cases from Mach 0.2 to 1.65:

| Case | Mach | Mesh size | Recorded GPU runtime | Status |
|---|---:|---:|---:|---|
| NACA 4412 LF1 / LF2 | 0.20 | 4,718 / 9,517 cells | 20.07 / 33.09 s total process time | $C_L$, $C_M$, and surface $C_p$ match XFOIL reference |
| RAE 2822 hybrid LF1 | 0.78 | 5,051 wedge/hex cells | 49,000 DP-LUR sweeps in 11.24 s | Solver verification passed |
| NASA CRM tetrahedral | 0.85 | 10,067,380 tetrahedra | 5,850 DP-LUR sweeps in 600.20 s (610.35 s total) | Throughput test; flow state not yet converged |
| NASA CRM Cartesian AMR | 0.85 | 21,523,036 active cells (from 31.48M octree cells) | 297 DP-LUR sweeps in 121.19 s (0.408 s/iter) | Large-grid throughput test; flow state not yet converged |
| 3D sphere regression | 1.65 | 7,680 polyhedral cells | 2.15 s total process time | Conservation and CPU/CUDA parity check |

On the 21.52M-cell CRM AMR grid, each DP-LUR iteration takes ~0.408 s on an RTX 4090. Building the full octree cut-cell mesh from the closed aircraft surface takes 8 min 52 s from scratch, or 5 min 20 s when refining from a cached parent grid.

## Mesh Representations

| Dim | Grid type | Application in this repo |
|---|---|---|
| 2D | Structured Cartesian & quadtree AMR cut cells | Multi-fidelity (LF1/LF2) airfoils with exact boundary segments and small-cell merging |
| 3D | Explicit owner/neighbor polyhedral faces | Unstructured tetrahedra, prism/hex hybrid slabs, and 6-block cubed-sphere grids |
| 3D | Procedural multi-block octree AMR cut cells | $8\times 8\times 8$ leaf blocks with exact surface clipping, cell agglomeration, and 2:1 refinement interfaces |

For the 3D octree AMR solver, regular interior block connectivity is evaluated procedurally on the fly; only block interfaces, physical boundaries, cut-face exceptions, and surface triangles are stored explicitly in device memory.

## Mesh and Flow Visualization

All surface and volume plots below are rendered directly on the true control-volume polygons, cut segments, or surface triangles without `scipy.interpolate.griddata` smoothing.

| 2D quadtree AMR cut cells | 3D cubed-sphere polyhedra |
|---|---|
| ![NACA4412 true AMR members, flow, and wall Cp](examples/2d_naca4412/showcase/naca4412_lf2.png) | ![Six-block sphere surface, exact plane cells, and mesh](examples/3d_sphere/showcase/sphere_m165_exact_mesh.png) |
| **Wedge/hex hybrid mesh** | **Body-fitted tetrahedral mesh** |
| ![RAE2822 wedge and hexahedron hybrid grid](examples/3d_rae2822_hybrid/showcase/rae2822_lf2_hybrid.png) | ![CRM exact tetrahedron-plane intersections](examples/3d_crm_tetra/showcase/crm_tetra_exact_section.png) |
| **3D Cartesian-octree AMR cut cells** | **Full 3D aircraft surface pressure** |
| ![CRM true leaf blocks and individual AMR cells](examples/3d_crm_amr/showcase/crm_amr_yz_x5_mesh.png) | ![CRM three-dimensional Cp perspectives](examples/3d_crm_amr/showcase/crm_amr_surface_cp_3d.png) |

Additional symmetry-plane, spanwise slice, and surface pressure plots are indexed in the [example showcase](examples/README.md#visual-showcase).

## Solver Paths

All solvers use a density-based finite-volume formulation (explicit RK3 or implicit DP-LUR) with low-Mach preconditioning for subsonic cases; pressure-based segregated schemes (SIMPLE/PISO) are not used.

| Path | Mesh topology | Discretization & time marching | Verification / validation |
|---|---|---|---|
| 2D airfoil | Cartesian cut cell (LF1/LF2) | Density-based, RK3, SA-neg, low-Mach preconditioning | Validated against XFOIL ($C_L, C_M, C_p$) on NACA 4412 |
| 3D sphere | 6-block cubed-sphere polyhedra | Explicit-face CPU/CUDA | Conservation, positivity, CPU/CUDA parity |
| 3D RAE 2822 | Mixed wedge/hex thin wing section | Explicit-face DP-LUR | LF1 run and CPU/CUDA parity |
| 3D CRM tetra | Body-fitted unstructured tetrahedra | Explicit-face DP-LUR | 10.07M-cell import, checkpoint restart, gradient & SA tests |
| 3D CRM AMR | Multi-block octree cut cell | Procedural DP-LUR | 7.66M, 14.01M, and 21.52M control-volume runs |

## Build and Quick Start

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

## Documentation

- [Theory and operation manual](docs/MANUAL.md)
- [Development log and history](docs/DEVELOPMENT_HISTORY.md)
- [Example matrix and runbooks](examples/README.md)
- [Mesh and flow visual showcase](examples/README.md#visual-showcase)
- [Current verification & convergence status](docs/STATUS.md)
- [Planned wall-reconstruction roadmap](docs/STATUS.md#next-development-phase)
- [Binary and mesh file formats](docs/FORMATS.md)
- [References and external datasets](docs/REFERENCES.md)

## Repository Layout

- `src/`, `pipeline/`, `quick_cfd.py` — 2D `quick_airfoil` Cartesian cut-cell workflow
- `solvers/` — [2D/3D solver executables and source map](solvers/README.md)
- `solvers/3d/` — Explicit-face `sphere3d` and procedural-AMR `amr3d` solvers, mesh importers, restart tools, and unit tests
- `tools/cartesian_amr/` — Octree builder, exact triangle-box clipping, cut-cell moments, cell agglomeration, and GPU table packing
- `tools/hybrid/` — 2D mixed tri/quad extrusion to 3D wedge/hex hybrid grids
- `examples/` — Reproducible runbooks and configuration files for each test case
- `validation/` — Post-processing and reference comparison scripts

## External Meshes and License

Raw NASA CRM and RAE 2822 archive grids are not stored in the repository; the scripts under `examples/` document the download URLs, checksums, and conversion commands. See [LICENSE_NOT_SELECTED.md](LICENSE_NOT_SELECTED.md) regarding current licensing status.
