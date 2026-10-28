# Mesh Data Structures and Reconstruction

This guide follows the current 3-D readers and residual operators. It explains
how node, cell, face, and AMR block data are related, and how that data reaches
the flux and reconstruction calculations. The 2-D airfoil path constructs its
Cartesian cut-cell mesh from an ordered boundary CSV; it does not read external
2-D volume meshes.

## 1. Two different 3-D mesh routes

| Route | Surface/volume input | Solver representation | Executable |
|---|---|---|---|
| Cartesian AMR cut-cell | Closed triangulated surface | Packed blocks, interfaces, cut faces, wall fragments and controls (`QCFDAM1`) | `amr3d` |
| Unstructured | Tetrahedra or mixed polyhedra converted to explicit faces | Cells, faces, adjacency and gradient terms (`QCFD3D1`) | `sphere3d` |

Both solve finite-volume equations, but their neighbor relationships are
represented differently. `sphere3d` stores an explicit face list and cell-face
adjacency. `amr3d` generates regular Cartesian block faces procedurally and
stores only inter-block and geometric exception faces.

## 2. Unstructured mesh: nodes, cells, and faces

### 2.1 Source connectivity versus solver connectivity

An imported volume mesh usually starts with **nodes** (coordinates) and
**elements** (tetrahedra, wedges, hexes, or other polyhedra). Element
connectivity answers which nodes form each cell. Two cells are neighbors when
they share a face, but the solver does not search all node lists during every
iteration. Preprocessing converts that topology into a face table and
cell-to-face adjacency.

In `QCFD3D1`, node IDs are generally absent from the flow-solver core. The
explicit-face mesh contains geometry sufficient for its finite-volume fluxes.
For tetrahedral node-based gradients, node coordinates and tetra connectivity
are retained separately in the `QCFDNC1` companion.

### 2.2 Face orientation and neighbor lookup

Each `Face` record has:

| Field | Meaning |
|---|---|
| `left` | Owner cell index; always present |
| `right` | Neighbor cell index for an interior face; negative for a boundary |
| `bc` | Boundary-condition code when `right < 0` |
| `x[3]` | Face centroid |
| `n[3]` | Unit normal directed outward from `left` toward `right` or the boundary |
| `area` | Face area |
| `distance` | Face-normal owner-neighbor distance, or owner-to-boundary distance |

For an interior face, the solver computes one numerical flux with this
orientation. It adds `+F A` to the owner residual and `-F A` to the neighbor
residual. This paired update is the basic local conservation mechanism. To
visit a cell's faces, the solver reads `Cell.first` and `Cell.count` into the
flat integer `ids` array. The corresponding face index is `ids[i]`; the face
record then gives the other cell directly through `left/right`.

```text
cell 12 (left/owner)  ---- face 87, normal n ---->  cell 29 (right/neighbour)
       residual += F*A                         residual -= F*A
```

Boundary faces use `right < 0`; `bc` selects the boundary treatment. A wall
face is also how surface pressure and viscous traction are integrated. A
reversed face must swap owner/neighbour and reverse its normal together; doing
only one of those changes gives the wrong flux sign.

### 2.3 QCFD3D1 record layout

The binary file starts with `MeshHeader` (`QCFD3D1` magic and counts), followed
by contiguous arrays in this order:

1. `Cell[n]`
2. `Face[nf]`
3. cell-face face IDs `int[ni]`
4. gradient `Term[nt]`

The current C++ records in `solvers/3d/mesh_format.h` are:

```cpp
struct Cell { double x[3], volume, dwall;
              int first, count, term_first, term_count; };
struct Face { int left, right, bc;
              double x[3], n[3], area, distance; };
struct Term { int source, face; double coefficient[3]; };
```

The text reader uses the same logical arrays. It reads the four counts first,
then cell rows, face rows, adjacency IDs, and gradient-term rows. The binary
reader detects the `QCFD3D1` magic and reads the packed arrays directly. Both
populate the same host arrays; the CUDA path then copies the required arrays to
the device. The exact reader is in `solvers/3d/solver.cu` near mesh loading.

### 2.4 Gradient and second-order face states

The default cell-centered gradient uses the preassembled `Term` records. For
cell `i`, each term names a source cell (or a boundary sentinel) and a vector
coefficient. In simplified form:

\[
\nabla q_i \approx \sum_{k\in stencil(i)} c_{ik}(q_k-q_i).
\]

When the node-gradient option is selected for supported tetrahedral inputs,
the `QCFDNC1` connectivity and `QCFDNG1` stencil data supply the node-based
weights and face-node information instead. It is not inferred from
`QCFD3D1`, because that file alone does not contain the source node topology.

For second-order mean-flow convection, the solver evaluates a limited
cell-centered gradient at the face centroid:

\[
q_{i,f}=q_i+\beta\,\phi_i\,\nabla q_i\cdot(x_f-x_i),
\]

where `x_i` is the cell centroid, `x_f` the face centroid, `beta` the selected
reconstruction blend, and `phi_i` a limiter in `[0,1]`. The limiter constrains
the face increment using neighboring values so the reconstruction does not
overshoot the local neighbor range. Interior faces reconstruct a state from
each side and pass both to HLLC. First-order mode uses the cell averages
directly. The code paths are `compute_gradient` and
`reconstruct_face_state` in `solvers/3d/solver.cu`.

The viscous face gradient has a separate correction: averaged cell gradients
are adjusted along the owner-neighbor displacement to match the observed
two-cell difference. The optional Mitchell node/face stencil is another
specialized route; it should not be conflated with the convective limiter.
SA convection remains first order in the retained general 3-D path unless an
explicit higher-order option is selected and verified.

## 3. AMR cut-cell mesh: topology and tables

### 3.1 From surface to control volumes

The AMR builder starts from a closed triangle surface. It creates an octree of
Cartesian leaf blocks (each currently 8x8x8 cells), classifies cells as fluid,
solid, or cut, computes cut-cell volume/centroid/open-face geometry, and
merges small fluid cells into larger agglomerated control volumes. Thus a
solver **control volume** can represent one ordinary Cartesian cell or a group
of procedural cells. Wall fragments retain the exact clipped solid boundary
used for no-slip fluxes and surface loads.

### 3.2 What the QCFDAM1 directory stores

`header.bin` contains magic, version, endianness, block size, root geometry,
and counts for the major arrays. The current reader then loads:

| File | Relationship/data |
|---|---|
| `blocks.bin` | Leaf-block origin, spacing, level and first procedural-cell index |
| `interfaces.bin` | Neighboring block pair, side and level difference (including 2:1 interfaces) |
| `boundaries.bin` | External block sides |
| `controls.bin` | Agglomerated control centroid, volume and wall distance |
| `procedural_control.bin` | Procedural Cartesian cell -> agglomerated control ID; negative means inactive/solid |
| `procedural_status.bin` | Per-procedural-cell geometry/status flags |
| `cut_faces.bin` | Open Cartesian face exceptions at cut or agglomerated geometry |
| `wall_faces.bin` | Exact wall-fragment centroid, normal, area and adjacent control |
| `wall_triangle.bin` | Wall fragment -> source surface triangle ID |
| `wall_distance.bin` | Wall distance for each control volume |
| `wall_counts.txt`, `wall_points.bin`, `wall_triangles.bin` | Original triangulated surface for loads and visualization |

Regular faces inside a block are not stored individually. Their locations,
areas, and axis normals follow from the block origin and spacing. The kernel
generates these faces, looks up each side through `procedural_control`, and
skips solid/inactive cells or internal faces whose two procedural cells map to
the same agglomerated control. At block boundaries it uses the interface
tables, including the level difference, to connect coarse and fine sides. Cut
faces and physical walls are read from their explicit exception tables.

### 3.3 AMR gradients and reconstruction

The AMR residual pass first accumulates face-area-weighted primitive values
for each control and divides by its fluid volume, giving a Green-Gauss
gradient. Procedural interior faces, inter-block faces, cut-face exceptions,
and boundary states contribute according to their geometry. The same-level
and 2:1 interface kernels use the mapped control IDs so flux is applied
conservatively to both participating controls.

With second-order reconstruction enabled, the state at a face is formed from
the control-center value and the stored gradient dotted with the center-to-face
offset. The increment is bounded by the neighboring value; density and
pressure candidates also have positivity/range guards. HLLC uses the two
reconstructed states on interior faces. At walls, no-slip and wall-pressure
handling follow the dedicated wall operator. This is second-order mean-flow
convection; it does not imply that every SA and viscous term is second order.

The AMR reader validates the `QCFDAM1` header, reads exact counts from it, then
loads the named arrays into host vectors. For CUDA runs it allocates device
arrays and copies the mesh tables once before iteration; the flow state and
residual work arrays then stay on the GPU. The implementation is in
`solvers/3d/amr_solver.cu`, where `load_exact` and `main` handle the reader,
`procedural_face` creates regular faces, and `reconstruct_face_state` evaluates
the face states.

## 4. Running the examples

```bash
# Self-contained 3-D sphere: generates an explicit-face mesh, then uses sphere3d.
bash examples/3d_sphere/run_smoke.sh build_cpu cpu

# 3-D CRM AMR: provide a closed surface, build QCFDAM1, then run amr3d.
export CRM_SURFACE=/path/to/CRM_full_closed.npz
bash examples/3d_crm_amr/build_meshes.sh build output/3d_crm_amr
bash examples/3d_crm_amr/run_upperwing.sh build cuda output/3d_crm_amr

# Unstructured CRM: provide/download CGNS, convert to QCFD3D1, then use sphere3d.
DOWNLOAD_NASA_CRM=1 bash examples/3d_crm_tetra/run.sh build cuda
```

Full steps and resource requirements are in the individual `RUNBOOK.md`
files. A mesh and its checkpoint must remain paired; `state.bin` contains
flow variables but no mesh identity.

## 5. Related format references

- [`FORMATS.md`](FORMATS.md): compact inventory of all binary formats.
- [`solvers/3d/mesh_format.h`](../solvers/3d/mesh_format.h): explicit-face,
  connectivity, gradient and agglomeration record declarations.
- [`solvers/3d/amr_mesh_format.h`](../solvers/3d/amr_mesh_format.h): AMR
  record declarations.
- [`examples/3d_crm_amr/RUNBOOK.md`](../examples/3d_crm_amr/RUNBOOK.md): full
  AMR construction, audits and solve process.
