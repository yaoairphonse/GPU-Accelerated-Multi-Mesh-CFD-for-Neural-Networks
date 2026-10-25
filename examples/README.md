# Example Matrix and Visual Showcase

Each case directory contains reproducible run scripts, solver configuration files, reference metadata, and unsmoothed mesh/flow plots.

| Example | Self-contained | Focus | Status |
|---|---|---|---|
| `3d_crm_amr` | requires CRM surface & GPU memory | 21.52M-cell procedural multi-block octree cut-cell aircraft | 5,000-iter run; wing $C_p$ RMSE 0.1003, $C_M$ mode active |
| `3d_crm_tetra` | downloads public NASA archive | 10.07M-cell unstructured body-fitted tetrahedral aircraft | Solver & restart verified; throughput test |
| `3d_sphere` | yes | 3D conservation, positivity, restart, CPU/CUDA parity ($M=1.65$) | Verification passed ($10^{-15}$ CPU/CUDA parity) |
| `3d_rae2822_hybrid` | requires archived NPZ inputs | Mixed wedge/hex transonic wing slab ($M=0.78$) | LF1 solver verification passed |
| `2d_naca4412` | yes | Multi-fidelity (LF1/LF2) low-Mach cut-cell airfoil ($M=0.20$) | Validated against XFOIL ($C_L, C_M, C_p$) |
| `2d_rae2822_amr` | geometry bundled | Transonic 2D Cartesian AMR cut-cell airfoil | Development test case |

Run all scripts from the repository root. Write solver outputs and checkpoints to `output/` or an external scratch directory rather than inside the source tree, and direct each restart continuation to a fresh output folder so input checkpoints are preserved.

## Visual Showcase

Every figure below is plotted directly from true mesh connectivity or procedural cell ownership without `griddata` smoothing.

| Case | Gallery |
|---|---|
| CRM Cartesian AMR (3D aircraft) | [3D aircraft $C_p$](3d_crm_amr/showcase/crm_amr_surface_cp_3d.png) · [surface $C_p$ & load history](3d_crm_amr/showcase/crm_surface_cp.png) · [wing $x$-$y$ Mach/$C_p$](3d_crm_amr/showcase/crm_xy_wing_mach_cp.png) · [full-wing $y$-$z$ Mach/$C_p$](3d_crm_amr/showcase/crm_yz_fullwing_mach_cp.png) · [symmetry $x$-$z$ Mach/$C_p$](3d_crm_amr/showcase/crm_xz_symmetry_mach_cp.png) · [AMR sections](3d_crm_amr/showcase/crm_amr_true_mesh_sections.png) · [wing cut-cell zoom](3d_crm_amr/showcase/crm_amr_yz_x5_mesh.png) · [notes](3d_crm_amr/showcase/README.md) |
| CRM tetrahedral (3D aircraft) | [exact wing section](3d_crm_tetra/showcase/crm_tetra_exact_section.png) · [notes](3d_crm_tetra/showcase/README.md) |
| 3D sphere ($M=1.65$) | [surface/field/mesh](3d_sphere/showcase/sphere_m165_exact_mesh.png) · [notes](3d_sphere/showcase/README.md) |
| RAE 2822 wedge/hex hybrid | [LF2 field/mesh/$C_p$](3d_rae2822_hybrid/showcase/rae2822_lf2_hybrid.png) · [notes](3d_rae2822_hybrid/showcase/README.md) |
| NACA 4412 Cartesian AMR (2D) | [LF1](2d_naca4412/showcase/naca4412_lf1.png) · [LF2](2d_naca4412/showcase/naca4412_lf2.png) · [notes](2d_naca4412/showcase/README.md) |
| RAE 2822 Cartesian AMR (2D) | [LF2 field/mesh/$C_p$](2d_rae2822_amr/showcase/rae2822_cartesian_amr_lf2.png) · [notes](2d_rae2822_amr/showcase/README.md) |
