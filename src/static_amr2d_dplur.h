#ifndef CFD_STATIC_AMR2D_DPLUR_H
#define CFD_STATIC_AMR2D_DPLUR_H

#include "static_amr2d.h"

namespace cfd {

// Steady local-pseudo-time controls.  This reference implementation is kept
// separate from physical-time SSP-RK3 so NNLCI production cases can switch
// integrators without changing their mesh, state, or output formats.
struct StaticAmrDplurOptions2D {
  long iterations = 100;
  int sweeps = 5;
  double cfl = 4;
  double omega = .8;
  double budget_seconds = 0; // zero disables the wall-clock stop
  bool block_coupling = true;
};

struct StaticAmrDplurReport2D {
  long accepted = 0;
  long rejected = 0;
  double seconds = 0;
  double final_cfl = 0;
  double residual[5] = {};
};

// Returns 0 after the requested updates, 2 on a clean budget stop, and 1 on
// invalid input or an unrecoverable nonphysical update.  Accepted updates are
// always retained; a rejected candidate never modifies state->u.
int static_amr2d_run_dplur(const StaticAmrMesh2D&, const Params2D&,
                           const StaticAmrDplurOptions2D&,
                           StaticAmrState2D*, StaticAmrDplurReport2D*,
                           std::string* error);

} // namespace cfd
#endif
