# AGENTS.md — This File Is for Your Agent Operation

This document is written for coding agents operating this repository. Follow
it before building, running, reviewing, or modifying the project.

## 1. Project scope

This repository contains several related finite-volume solver paths rather
than one universal executable:

| Path | Dimension | Mesh | Entry point |
|---|---:|---|---|
| Fast airfoil workflow | 2-D | Cartesian/quadtree AMR cut cells | `simple_airfoil.py` -> `quick_airfoil` |
| Full airfoil workflow | 2-D | Cartesian/quadtree AMR cut cells | `quick_cfd.py` -> `quick_airfoil` |
| Procedural aircraft AMR | 3-D | Cartesian octree AMR cut cells | `amr3d` |
| General unstructured solver | 3-D | tetrahedral, wedge/hex hybrid, or polyhedral | `sphere3d` |

Do not look for a single `--dimension 2d|3d` switch. Select the entry point
from the table above. For a new user asking to run an airfoil, use the native
2-D path; do not extrude it into a thin 3-D mesh.

## 2. Required reading

Read these files before taking action:

1. `README.md` for the project status and solver matrix;
2. `docs/SIMPLE_AIRFOIL_INTERFACE.md` for the locked 2-D interface;
3. the matching `examples/<case>/RUNBOOK.md` before running a case;
4. `docs/STATUS.md` before making accuracy or validation claims.

Use `examples/README.md` to determine whether a case is self-contained. A
plotted showcase is evidence of an earlier run, not proof that its large input
mesh or checkpoint is stored in Git.

## 3. Default first operation: self-contained 2-D airfoil

Unless the user requests another case, verify the repository with the bundled
NACA 4412 LF1 CPU DP-LUR case. It requires no external mesh and avoids a CUDA
dependency.

Create the environment and build from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

cmake -S . -B build_cpu \
  -DCMAKE_BUILD_TYPE=Release \
  -DCFD_CUDA=OFF \
  -DCFD_OPENMP=ON
cmake --build build_cpu -j 8 --target quick_airfoil
ctest --test-dir build_cpu --output-on-failure
```

Then run:

```bash
python3 simple_airfoil.py \
  --geometry examples/2d_naca4412/geometry.csv \
  --algorithm dplur \
  --mesh-cells 4700 \
  --mach 0.2 \
  --aoa 1.5 \
  --pressure-pa 101325 \
  --temperature-k 288.15 \
  --out output/naca4412_lf1
```

`simple_airfoil.py` does not compile the executable. Build first, or set
`QUICK_CFD_EXECUTABLE` to an already-built `quick_airfoil` binary. Always use
a new output directory rather than overwriting a prior run.

The primary machine-readable answer is:

```text
output/naca4412_lf1/simple_summary.json
```

Also inspect `pipeline.json`, the per-level `result.json`, convergence history,
surface pressure, and the generated HTML index before reporting success.

## 4. Interpreting the simple 2-D controls

- `--algorithm dplur` selects the CPU matrix-free DP-LUR reference.
- `--algorithm rk3` selects the CUDA RK3 path and requires a CUDA build.
- `--mesh-cells 4700` selects LF1; the NACA 4412 case has about 4,718 actual
  control volumes.
- `--mesh-cells 9500` selects LF2; the NACA 4412 case has about 9,517 actual
  control volumes.
- Geometry must be CSV with header `x,y` or ordered `x,y,patch`.
- Pressure and temperature use SI units; angle of attack is in degrees.

The interface deliberately locks spatial order, turbulence model,
preconditioning, CFL, relaxation, and convergence gates. Do not tune those
settings when the user's goal is merely to run the standard case.

Reference workstation timings are approximately 7.7 seconds for LF1 CPU
DP-LUR and 48.3 seconds for LF2 CPU DP-LUR. Timings are hardware dependent.
For the retained NACA 4412 LF1 DP-LUR run, the expected neighborhood is about
4,718 controls, 1,900 pseudo-iterations, and
`CL/CD/CM = 0.6600/0.0068/+0.1073`. Investigate material deviations; do not
hard-code these values as pass/fail results on a modified geometry.

Lift, pitching moment, and surface-pressure trends are the supported 2-D
accuracy outputs. Drag is exported but the rapid cut-cell viscous drag is not
yet an aerodynamic validation claim.

## 5. CUDA airfoil operation

Use CUDA only when the user requests it or the device is available. Check
`nvidia-smi` first and select the architecture appropriate to that GPU; do not
blindly assume RTX 4090 architecture 89 on another machine.

Example for an RTX 4090:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCFD_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=89
cmake --build build -j 8 --target quick_airfoil
ctest --test-dir build --output-on-failure
```

Run the same simple command with `--algorithm rk3`. Do not silently fall back
from a requested CUDA run to CPU; report the unavailable device or build
dependency clearly.

## 6. Three-dimensional case routing

- `examples/3d_sphere` is the self-contained 3-D regression case. Its mesh is
  generated locally. Prefer it for an initial 3-D smoke test.
- `examples/3d_crm_tetra` requires the public NASA CGNS archive. Download is
  opt-in through `DOWNLOAD_NASA_CRM=1` or an explicit `CRM_CGNS` path.
- `examples/3d_crm_amr` requires a closed CRM surface, several gigabytes of
  generated data, and a high-memory GPU for the retained 21.52M-control case.
- `examples/3d_rae2822_hybrid` requires external archived NPZ topology and
  true-mesh inputs that are not redistributed in this repository.

The retained full CRM AMR mesh and 11,900-step checkpoint are not in GitHub.
Do not claim that cloning this repository alone reproduces that exact restart.
Do not download, generate, copy, or commit multi-gigabyte meshes/checkpoints
without explicit user direction.

## 7. Numerical and reporting guardrails

- Distinguish verification, benchmarking, and aerodynamic validation.
- XFOIL is a numerical comparison, not experimental validation.
- The CRM steady pseudo-time trajectory is not fixed-point converged. Its
  observed low-to-low cycle mean differs from NASA USM3D by approximately
  `-3.39%/+4.25%/+15.51%` in `CD/CL/CM`; individual iteration snapshots can
  be closer or farther away.
- A pseudo-iteration average is not a physical-time URANS average.
- Residual smoothing or prolongated corrections are not full FAS multigrid.
- Preserve the documented moment reference point and sign convention when
  comparing `CM`.
- Never report a successful run from coefficients alone. Check positivity,
  residual behavior, convergence status, rejected updates, and load history.

## 8. Repository operation

- Keep generated builds and outputs in `build*` and `output/`; do not add them
  to source control.
- Preserve user changes and inspect `git status` before editing.
- Do not change locked example inputs merely to make a test pass.
- For code changes, build the affected CPU/CUDA target and run the relevant
  CTest plus at least one small case.
- Keep published documentation, source comments, configuration keys, and CLI
  text in English.
- Do not commit or push unless the user explicitly asks.
- Do not make quantitative claims stronger than the evidence recorded in the
  runbooks and `docs/STATUS.md`.

When handing work back, state the solver path, mesh/control count, algorithm,
backend, iteration count, wall-clock time, convergence status, output path,
and any missing external input.
