#ifndef CFD_STATIC_AMR2D_DIAGNOSTICS_H
#define CFD_STATIC_AMR2D_DIAGNOSTICS_H

#include "static_amr2d.h"

namespace cfd {
struct StaticAmrAuditOptions2D {
  double angle_radians = 0;
  double dynamic_pressure = 1;
  double reference_length = 1;
  double moment_x = 0;
  double moment_y = 0;
  int bins = 24;
};

// Writes PREFIX.face_budget.csv and PREFIX.wall_load_budget.csv.  The same
// face flux, wall pressure, viscous stress, and SA operators used by the
// solver are reused; no scattered interpolation or surrogate force is used.
int static_amr2d_write_budgets(const StaticAmrMesh2D&,const Params2D&,
                               const StaticAmrState2D&,
                               const StaticAmrAuditOptions2D&,
                               const std::string& prefix,std::string* error);
}
#endif
