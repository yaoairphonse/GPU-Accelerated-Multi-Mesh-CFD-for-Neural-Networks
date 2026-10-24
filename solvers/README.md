# Solver map

The repository contains three production/development solver entry points.  The
native 2-D implementation predates the `solvers/` directory and remains under
the root `src/` tree to preserve its build and Python-pipeline interface.

| Executable | Dimension | Mesh representation | Source entry point |
|---|---:|---|---|
| `quick_airfoil` | 2-D | structured/quadtree Cartesian AMR cut cells | [`src/static_amr2d*`](../src/) and [`solvers/2d/README.md`](2d/README.md) |
| `sphere3d` | 3-D | explicit owner/neighbour polyhedral faces, including tetra/wedge/hexa | [`solvers/3d/solver.cu`](3d/solver.cu) |
| `amr3d` | 3-D | procedural multi-block Cartesian-octree AMR cut cells | [`solvers/3d/amr_solver.cu`](3d/amr_solver.cu) |

The explicit-face 3-D solver can also consume a one-cell-thick extrusion of a
2-D mesh with symmetry caps.  The RAE2822 wedge/hexa example demonstrates that
quasi-2-D path.  For a genuinely two-dimensional calculation, `quick_airfoil`
is the smaller and faster implementation.
