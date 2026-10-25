# NASA CRM Cartesian-AMR Showcase

Flow slices below map `QCFDAM1` control-volume states directly through the procedural block and cell ownership tables, with solid boundaries rendered from exact plane intersections against the CRM surface triangles (no scattered `griddata` interpolation).

- [Full-aircraft 3D surface $C_p$ perspectives](crm_amr_surface_cp_3d.png) — conservative cut-wall pressure mapped onto all 194,980 surface triangles
- [Multi-view surface $C_p$ projections and load history](crm_surface_cp.png)
- [Horizontal wing slice ($x$-$y$ plane at $z/\text{MAC}=0.80$, $C_p$ and Mach)](crm_xy_wing_mach_cp.png)
- [Full-span crossflow slice ($y$-$z$ plane at $x/\text{MAC}=5.0$, $C_p$ and Mach)](crm_yz_fullwing_mach_cp.png)
- [Near-symmetry slice ($x$-$z$ plane at $y \to 0^+$, $C_p$ and Mach)](crm_xz_symmetry_mach_cp.png)
- [Octree-AMR block and cut-cell cross-sections](crm_amr_true_mesh_sections.png)
- [Spanwise leaf-block overview and wing cut-cell zoom at $x/\text{MAC}=5.0$](crm_amr_yz_x5_mesh.png)

The grid contains 21,523,036 active control volumes with a finest Cartesian spacing of $0.00390625\,\text{MAC}$. Figures reflect cumulative iteration 5,000 (4,000–5,000 window average: $C_L = 0.5237$, $C_D = 0.0263$, $C_M = -0.0245$, wing $C_p$ RMSE $0.1003$).
