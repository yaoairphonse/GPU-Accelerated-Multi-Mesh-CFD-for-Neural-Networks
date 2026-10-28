# NASA CRM Cartesian-AMR Showcase

Flow slices below map `QCFDAM1` control-volume states directly through the procedural block and cell ownership tables, with solid boundaries rendered from exact plane intersections against the CRM surface triangles (no scattered `griddata` interpolation).

- [Full-aircraft 3D surface $C_p$ perspectives](crm_amr_surface_cp_3d.png) — conservative cut-wall pressure mapped onto all 194,980 surface triangles
- [Multi-view surface $C_p$ projections and load history](crm_surface_cp.png)
- [Horizontal wing slice ($x$-$y$ plane at $z/\text{MAC}=0.80$, $C_p$ and Mach)](crm_xy_wing_mach_cp.png)
- [Full-span crossflow slice ($y$-$z$ plane at $x/\text{MAC}=5.0$, $C_p$ and Mach)](crm_yz_fullwing_mach_cp.png)
- [Near-symmetry slice ($x$-$z$ plane at $y \to 0^+$, $C_p$ and Mach)](crm_xz_symmetry_mach_cp.png)
- [Cumulative load and residual history through 11,900 pseudo-iterations](crm_amr_convergence_11900.png)
- [Octree-AMR block and cut-cell cross-sections](crm_amr_true_mesh_sections.png)
- [Spanwise leaf-block overview and wing cut-cell zoom at $x/\text{MAC}=5.0$](crm_amr_yz_x5_mesh.png)

The grid contains 21,523,036 active control volumes with a finest Cartesian spacing of $0.00390625\,\text{MAC}$. Refreshed flow and pressure figures use cumulative pseudo-iteration 11,900: $C_D/C_L/C_M = 0.029967/0.572087/-0.041837$. The three volume slices contain zero uncovered domain or unfilled fluid pixels at the recorded raster resolution.

This remains a development state rather than a converged CRM load result. The full observed low-to-low trajectory cycle averages $C_D/C_L/C_M = 0.028936/0.565420/-0.052798$ versus NASA USM3D $0.02995/0.54239/-0.06249$. It is a steady pseudo-iteration diagnostic average, not a physical-time URANS average.
