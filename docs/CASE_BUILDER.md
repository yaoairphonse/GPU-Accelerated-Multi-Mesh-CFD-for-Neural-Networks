# Offline case and input builder

Open [`tools/cfd_case_builder.html`](../tools/cfd_case_builder.html) in a web
browser after cloning the repository.  It is a self-contained offline page:
it has no server, package, CDN, telemetry, or network dependency.  It creates
inputs and commands but does not execute CFD in the browser.

## Supported routes

| Route in the page | Solver | Mesh input | Native configuration |
|---|---|---|---|
| 2-D native Cartesian/quadtree AMR airfoil | `quick_cfd.py` -> `quick_airfoil` | airfoil CSV | downloadable partial `case.json` plus command |
| 3-D native procedural octree-AMR cut cell | `amr3d` | packed `QCFDAM1` directory | explicit CLI plus manifest |
| 3-D external tetra/hybrid/polyhedral | `sphere3d` | prepared `QCFD3D1`/explicit-face file | explicit CLI plus manifest |

The page deliberately does not offer unsupported combinations. In particular,
the current native 2-D path does not read an external volume mesh, the 2-D
driver does not expose an inviscid Euler switch, and the 3-D production paths
use the six-variable SA-neg mechanical-viscous model. The native 2-D DP-LUR
integrator is available on CPU and CUDA; CUDA uses persistent device state and
CUDA Graph batches.

## Inputs represented

The form records:

- solver dimension and mesh topology;
- native AMR fidelity or prepared external mesh path;
- turbulence/viscosity options that exist on the selected path;
- RK3 explicit or DP-LUR implicit marching;
- spatial order, flux/preconditioning where available, CFL, DP-LUR sweeps,
  relaxation, and block coupling;
- iteration, convergence/stop-time, report interval, and wall-clock limits;
- uniform-freestream or 3-D restart initialization;
- Mach, angle of attack, static pressure/temperature for 2-D, and 3-D
  Reynolds/mu-star/reference-load conventions;
- build commands, run command, required preprocessing, output expectations,
  and known limitations.

An important interface distinction is preserved: 2-D RK3 stops through the
configured convergence gate, nondimensional stop time, or wall-clock budget;
it does not currently accept a direct maximum-iteration setting.  The 2-D
DP-LUR and both 3-D executables expose an iteration limit.

For the retained Mach 0.2 NACA 4412 DP-LUR setup, CFL 0.5 is the verified
portable LF1/LF2 selection. CFL 4 passed LF1 but produced a non-decaying LF2
residual cycle even though its loads changed only slightly; the page therefore
warns against treating that high-CFL endpoint as steady.

## Downloads

The page can download three products:

1. `<case>.manifest.json` — the normalized, agent-readable handoff;
2. `<case>.case.json` — a native partial configuration for the 2-D
   `quick_cfd.py --config` path (2-D only);
3. `run_<case>.sh` — build, test, and run commands with shell-quoted paths.

The manifest declares `quick_cfd_case_manifest_v1` and follows
[`tools/case_manifest.schema.json`](../tools/case_manifest.schema.json).  It
contains both an argument array and a display shell command.  Agents should
prefer the argument array when programmatically launching a process and
should still verify paths, external files, available GPU architecture, and
physical appropriateness before execution.

## Agent procedure

When a user supplies a generated manifest:

1. check `schema_version`, `case.solver_path`, and `case.dimension`;
2. read the corresponding example runbook;
3. verify every local path and external prerequisite without changing it;
4. inspect the `limitations` array;
5. when `native_configuration` is not null, write it verbatim to the case
   filename referenced after `--config` in `run.argv`;
6. build and run tests using `build.commands`;
7. use a new output directory;
8. execute `run.argv` without reinterpreting shell quoting;
9. report mesh/control count, backend, integrator, iterations, time,
   convergence state, loads, output path, and any warning.

The executable and Python configuration validator are authoritative.  A
manifest records an intended operation; it is not proof of a valid mesh,
stable numerical setup, convergence, or aerodynamic accuracy.

## Large meshes and external inputs

The 3-D AMR selection expects an already packed mesh directory.  Building an
octree from a closed surface remains the separate process documented in
`examples/3d_crm_amr/RUNBOOK.md`.  The external unstructured selection expects
a prepared explicit-face solver mesh; raw NASA CGNS and archived RAE NPZ files
must first pass their case-specific converters and audits.  Do not treat a
browser-selected filename as converted or validated data.
