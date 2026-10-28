# Installation and Operation Guide

This repository contains several solver entry points.  It is **not** one
universal executable with a `--dimension` switch.  Select the route that
matches the requested geometry and mesh:

| Goal | Dimension | Entry point | First example |
|---|---:|---|---|
| Fast airfoil/NNLCI input | 2-D | `simple_airfoil.py` -> `quick_airfoil` | `examples/2d_naca4412` |
| Fully configurable airfoil | 2-D | `quick_cfd.py` -> `quick_airfoil` | `examples/2d_naca4412` |
| Self-contained 3-D regression | 3-D | `sphere3d` | `examples/3d_sphere` |
| Cartesian aircraft AMR | 3-D | `amr3d` | `examples/3d_crm_amr` |
| Imported tetrahedral/polyhedral mesh | 3-D | `sphere3d` | `examples/3d_crm_tetra` |
| Wedge/hex hybrid mesh | 3-D | `sphere3d` | `examples/3d_rae2822_hybrid` |

The repository includes source code, runbooks, showcase figures, and compact
reference metadata.  It does not include the multi-gigabyte CRM packed mesh,
NASA archive files, RAE2822 NPZ topology, or restart checkpoints.

## 1. Clone and install dependencies

Run these commands from a Linux or WSL terminal:

```bash
git clone https://github.com/yaoairphonse/GPU-Accelerated-Multi-Mesh-CFD-for-Neural-Networks.git
cd GPU-Accelerated-Multi-Mesh-CFD-for-Neural-Networks

python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The Python requirements are used by the mesh builders, plotting tools, and
2-D pipeline.  A C++17 compiler and CMake 3.18 or newer are required.  CUDA
is optional for CPU operation and required for the CUDA backend.

## 2. Build the CPU version

This builds the 2-D solver and the 3-D `sphere3d`/`amr3d` tools without CUDA:

```bash
cmake -S . -B build_cpu \
  -DCMAKE_BUILD_TYPE=Release \
  -DCFD_CUDA=OFF \
  -DCFD_OPENMP=ON
cmake --build build_cpu --parallel 8
ctest --test-dir build_cpu --output-on-failure
```

The first smoke test should be the bundled 2-D CPU DP-LUR case:

```bash
python3 simple_airfoil.py \
  --geometry examples/2d_naca4412/geometry.csv \
  --algorithm dplur \
  --mesh-cells 4700 \
  --mach 0.2 --aoa 1.5 \
  --pressure-pa 101325 --temperature-k 288.15 \
  --out output/naca4412_lf1_cpu
```

The main result is `output/naca4412_lf1_cpu/simple_summary.json`.  The same
directory contains the surface pressure, loads, convergence history, and an
HTML index.  Always use a new output directory for a new run.

## 3. Build and run CUDA

Check the GPU before configuring CMake:

```bash
nvidia-smi
```

Then build with the architecture matching the installed GPU.  `89` is the
RTX 4090 setting used for the published timing measurements; it is not a
universal default:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCFD_CUDA=ON \
  -DCFD_OPENMP=OFF \
  -DCMAKE_CUDA_ARCHITECTURES=89
cmake --build build --parallel 8
ctest --test-dir build --output-on-failure
```

For the minimal interface, `--algorithm rk3` selects the CUDA RK3 path:

```bash
python3 simple_airfoil.py \
  --geometry examples/2d_naca4412/geometry.csv \
  --algorithm rk3 \
  --mesh-cells 4700 \
  --mach 0.2 --aoa 1.5 \
  --pressure-pa 101325 --temperature-k 288.15 \
  --out output/naca4412_lf1_cuda
```

For the verified persistent CUDA block-DP-LUR route, use the explicit
configuration rather than the minimal six-input wrapper:

```bash
python3 quick_cfd.py run \
  --points examples/NACA4412_polygon.csv \
  --point-format polygon \
  --config config/naca4412_m02_order2_dplur_cuda.json \
  --levels LF1 LF2 \
  --backend cuda \
  --executable build/quick_airfoil \
  --out output/naca4412_dplur_cuda
```

Open `output/naca4412_dplur_cuda/index.html` after completion.  The retained
RTX 4090 measurements are approximately 3.22 s for LF1 (4,718 controls) and
6.25 s for LF2 (9,517 controls).  These are hardware- and case-specific
measurements, not a timing guarantee.

## 4. Run the self-contained 3-D case

The sphere case generates its own mesh and is the correct first 3-D test:

```bash
bash examples/3d_sphere/run_smoke.sh build_cpu cpu
```

For a CUDA build:

```bash
bash examples/3d_sphere/run_smoke.sh build cuda
```

The result is written below `output/3d_sphere_cpu/` or
`output/3d_sphere_cuda/`.  This is a numerical verification case, not a
sphere-drag validation benchmark.

## 5. Run the 3-D CRM AMR path

The CRM AMR path uses only a closed aircraft surface and builds a new
Cartesian octree/cut-cell mesh.  It does not use the NASA volume mesh.  Set
`CRM_SURFACE` to a closed CRM NPZ/STL surface, then follow the dedicated
runbook:

```bash
export CRM_SURFACE=/path/to/CRM_full_closed.npz
bash examples/3d_crm_amr/build_meshes.sh build output/3d_crm_amr
bash examples/3d_crm_amr/run_upperwing.sh build cuda output/3d_crm_amr
```

This is a multi-gigabyte workflow.  Read
[`examples/3d_crm_amr/RUNBOOK.md`](../examples/3d_crm_amr/RUNBOOK.md) before
starting it.  The retained 21.52M-control CRM checkpoint is not in GitHub;
cloning the repository alone does not reproduce that restart state.

## 6. Run imported tetrahedral or hybrid meshes

The tetrahedral CRM example downloads the public NASA archive only when
explicitly requested:

```bash
DOWNLOAD_NASA_CRM=1 \
bash examples/3d_crm_tetra/run.sh build cuda
```

Alternatively set `CRM_CGNS=/path/to/file.cgns`.  The imported mesh is
converted locally before `sphere3d` runs.

The RAE2822 hybrid path requires archived NPZ topology and true-mesh inputs:

```bash
export RAE_TOPOLOGY=/path/to/LF1_full_topology.npz
export RAE_TRUE_MESH=/path/to/LF1_true_mesh.npz
bash examples/3d_rae2822_hybrid/run.sh build cuda
```

Those files are intentionally not redistributed in this repository.

## 7. How to interpret outputs

Every run writes a result JSON and a provenance/configuration record.  Before
reporting a result, inspect:

1. `result.json` for exit status, rejected updates, and convergence state;
2. `history.csv` for residual and load trends;
3. `output/loads.json` or the exported aerodynamic table for force/moment
   components;
4. the generated HTML or showcase plot for the actual mesh/field geometry.

The current solver is density-based finite volume with HLLC, RK3 or DP-LUR,
and low-Mach preconditioning.  A pressure-based SIMPLE/PISO solver is **not**
included yet.  The CRM results are pseudo-time development diagnostics and
must not be described as converged NASA validation data.

## 8. Browser case builder

For a human- or agent-readable setup form, open
[`tools/cfd_case_builder.html`](../tools/cfd_case_builder.html) locally in a
browser.  It generates a case manifest and run script; it does not execute
CFD in the browser.  Read the selected example runbook and verify all external
paths before executing the generated command.

## 9. AMR and unstructured mesh formats

The 3-D paths use two different mesh representations. Select one based on the
available input and desired workflow:

| Route | What you provide | What the solver reads | Main preparation tool |
|---|---|---|---|
| Procedural Cartesian AMR (`amr3d`) | Closed, manifold aircraft surface (`.npz` or `.stl`) | A packed `QCFDAM1` directory | `tools/cartesian_amr/` pipeline |
| Unstructured (`sphere3d`) | Tetra, wedge/hex, or polyhedral mesh converted to project format | `QCFD3D1` text or binary explicit-face mesh | `prepare_tetra` or `tools/hybrid/convert_rae2822_npz.py` |

### AMR: how to use it

AMR starts from a watertight surface, builds an isotropic Cartesian octree,
clips cells against the surface to create cut cells, merges small fluid cells,
and packages the solver tables. The CRM scripts show the complete workflow:

```bash
export CRM_SURFACE=/path/to/CRM_full_closed.npz
bash examples/3d_crm_amr/build_meshes.sh build output/3d_crm_amr
bash examples/3d_crm_amr/run_upperwing.sh build cuda output/3d_crm_amr
```

`build_meshes.sh` writes intermediate planning and audit data under the case
directory and produces `output/3d_crm_amr/upperwing_l10/packed/`. The solver
then receives that directory with `--mesh`. The full input requirements,
memory expectations, audit commands, and continuation procedure are in the
[CRM AMR runbook](../examples/3d_crm_amr/RUNBOOK.md). A surface file alone is
not a solver-ready AMR mesh.

### AMR: what the packed format contains

`QCFDAM1` is a directory of versioned binary tables, not a single mesh file.
`header.bin` contains magic/version/endian tags and counts for blocks,
procedural cells, agglomerated controls, interfaces, boundaries, cut faces,
and wall faces. The reader uses those counts to load the corresponding
`blocks.bin`, `interfaces.bin`, `boundaries.bin`, `controls.bin`,
`procedural_control.bin`, `procedural_status.bin`, `cut_faces.bin`,
`wall_faces.bin`, `wall_triangle.bin`, and `wall_distance.bin` tables. Surface
points and triangles are stored in `wall_points.bin` and
`wall_triangles.bin`, with counts in `wall_counts.txt`.

The `8x8x8` leaf blocks keep regular interior connectivity procedural; the
tables store interfaces and geometric exceptions explicitly. Do not rename,
omit, hand-edit, or mix files from different packed mesh builds. See
[`solvers/3d/amr_mesh_format.h`](../solvers/3d/amr_mesh_format.h) for the
record layouts and the runbook for the packing stages.

### Unstructured: how to use it

`sphere3d` reads explicit owner/neighbour faces. For the NASA tetrahedral CRM,
the run script obtains or accepts a CGNS source, converts it to project mesh
tables, and then runs the solver:

```bash
DOWNLOAD_NASA_CRM=1 bash examples/3d_crm_tetra/run.sh build cuda
```

You can instead set `CRM_CGNS=/path/to/file.cgns`. For the RAE2822 wedge/hex
example, provide both archived NPZ inputs:

```bash
export RAE_TOPOLOGY=/path/to/LF1_full_topology.npz
export RAE_TRUE_MESH=/path/to/LF1_true_mesh.npz
bash examples/3d_rae2822_hybrid/run.sh build cuda
```

The converted `QCFD3D1` mesh stores cell centroids, volumes and wall
distances, face ownership/neighbour or boundary codes, face centroids,
oriented area normals, areas and distances, plus cell-face and gradient-term
tables. `sphere3d` checks the first eight bytes: `QCFD3D1` selects the binary
reader; otherwise it reads the documented whitespace-delimited text form.
Both readers populate the same in-memory cell, face, adjacency and gradient
arrays before CPU or CUDA execution. Tetrahedral connectivity and boundary
triangles needed by node-based gradients are separate `QCFDNC1` connectivity
and `QCFDNG1` stencil companion files; they are not embedded in `QCFD3D1`.

The solver does not directly parse arbitrary CGNS, Gmsh, Fluent, SU2, or
VTK files. Convert source data into the supported representation first. The
format index and binary record declarations are in
[`docs/FORMATS.md`](FORMATS.md) and
[`solvers/3d/mesh_format.h`](../solvers/3d/mesh_format.h).
