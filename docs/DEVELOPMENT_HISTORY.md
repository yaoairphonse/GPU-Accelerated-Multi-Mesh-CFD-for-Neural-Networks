# Step-by-step development history

1. **Regular-grid Euler core.** Established conservative 1-D/2-D Euler
   states, Riemann fluxes, limiters, and explicit time integration.
2. **Immersed geometry.** Added signed-distance/immersed-boundary experiments,
   then replaced qualitative wall masking with conservative polygon cut cells.
3. **LF1/LF2 static AMR.** Added two fast Cartesian fidelity levels, exact
   apertures and wall segments, small-cell stabilization, integrated loads, and
   structured exports for NNLCI.
4. **Viscous RANS.** Added compressible mechanical viscous terms, Sutherland
   or constant viscosity, wall distance, and SA-neg on CPU/CUDA.
5. **Low-Mach work.** Added density-based Turkel-style preconditioning and
   matched CPU/CUDA operators. This did not introduce a pressure-based solver.
6. **Second-order 2-D accuracy.** Added agglomerated-control least-squares
   primitive gradients, multidimensional limiting, reconstructed face/wall
   pressure, and the CL/CM/Cp convergence gate. NACA4412 became the retained
   fast accuracy case.
7. **True three-dimensional solver.** Generalized the residual to explicit
   owner/neighbour faces and six conserved variables. Built a six-block
   cubed-sphere case to verify closed 3-D conservation, the third velocity
   component, positivity, restart, and CPU/CUDA agreement.
8. **NASA CRM tetra import.** Imported the public 10.07M-cell CGNS grid,
   preserved ordering and zones, built binary solver/connectivity tables, wall
   distance, restart, surface loads, and diagnostics.
9. **Convergence and accuracy tools.** Added limited reconstruction,
   matrix-free block DP-LUR, adaptive CFL, persistent-reversal detection,
   metric-best checkpoints, rollback, node gradients, Mitchell face gradients,
   SA residual audits, and an agglomeration residual-smoothing prototype.
10. **Mixed hybrid consumer.** Converted RAE2822 triangles/quads into wedges
    and hexahedra without tetrahedral splitting, proving that the explicit-face
    solver consumes mixed 3-D cells.
11. **Direct CRM Cartesian AMR.** Extracted the closed CRM surface, created
    balanced 8-cell octree blocks, exact wall clipping, cut volumes/apertures,
    small-cell agglomeration, procedural regular residuals, conservative 2:1
    subfaces, wall distance, and the compact QCFDAM1 format.
12. **AMR accuracy refinement.** Progressed from 7.66M controls to 14.01M
    all-wall refinement and 21.52M targeted upper-wing Level-10 refinement.
    Geometry, uniform-flow, affine-gradient, and CPU/CUDA audits passed.
13. **Warm fine-grid initialization.** Added exact volume-weighted state
    prolongation between nested AMR meshes. The 14.0M-to-21.5M transfer
    preserved global conserved integrals to roundoff and avoided a full cold
    start.
14. **Honest result packaging.** Added direct leaf-cell slices, true square
    mesh sections, original-triangle surface Cp, canonical runbooks, status
    labels, and non-destructive cache/release manifests.
15. **Pressure-only wall probe isolation.** Added a fixed-donor, bounded wall-
    normal pressure extrapolation behind an opt-in switch.  Uniform-flow and
    CPU/CUDA tests passed, but a same-checkpoint CRM A/B changed wing Cp RMSE by
    only 0.018%; it was rejected as the primary accuracy correction.
16. **Bounded SA face backbone.** Added neighbor-bounded SA convection,
    positive bounded nonorthogonal diffusion, and a consistent Rusanov/two-
    point DP-LUR action.  It lowered the 300-update SA residual about 1.77x but
    left the load trajectory unchanged, so it remains opt-in.
17. **Procedural mean-flow block DP-LUR.** Added and finite-difference-tested a
    matrix-free five-variable Euler normal-flux Jacobian action.  CPU/CUDA one-
    update agreement was 1.87e-16.  On 21.52M controls it improved residuals
    only 1--2% at 10.4% additional compute cost and did not remove load drift.
18. **5,000-iteration CRM assessment.** Continued the 1,278-update checkpoint
    by 3,722 updates.  Global mean/SA residuals reached 6.99e-5/3.28e-4 and
    wing Cp RMSE reached 0.10034.  The iteration-weighted CL/CD/CM average over
    cumulative iterations 4,000--5,000 is
    0.523661/0.026297/-0.024508, but the last-500 average differs materially;
    no steady or physical-time convergence is claimed.

The remaining numerical milestone is to localize the low-frequency CRM mode
by face family and shock/sectional-moment tracking, then distinguish pseudo-
time behavior from physical unsteadiness before a targeted 30--40M grid study.
