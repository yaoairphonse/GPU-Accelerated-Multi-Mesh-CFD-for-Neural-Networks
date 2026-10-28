#ifndef CFD_STATIC_AMR2D_RECONSTRUCTION_H
#define CFD_STATIC_AMR2D_RECONSTRUCTION_H

#include "static_amr2d.h"

namespace cfd {
#ifdef __CUDACC__
#define AMR_RECON_HD __host__ __device__
#else
#define AMR_RECON_HD
#endif
AMR_RECON_HD inline double amr_face_x(const AmrFace2D& f) { return f.axis==0?f.coordinate:f.centroid; }
AMR_RECON_HD inline double amr_face_y(const AmrFace2D& f) { return f.axis==1?f.coordinate:f.centroid; }
AMR_RECON_HD inline double amr_wall_x(const WallFragment2D& w) { return .5*(w.x0+w.x1); }
AMR_RECON_HD inline double amr_wall_y(const WallFragment2D& w) { return .5*(w.y0+w.y1); }
AMR_RECON_HD inline double amr_point_value(const double* primitive,const double* slopes,size_t n,
                                            int variable,size_t cell,const AmrCell2D& c,double x,double y) {
  return primitive[variable*n+cell]+slopes[(2*variable)*n+cell]*(x-c.cx)+slopes[(2*variable+1)*n+cell]*(y-c.cy);
}
AMR_RECON_HD inline double amr_limiter(double center,double lower,double upper,double delta) {
  if(delta>0&&center+delta>upper)return (upper-center)/delta;
  if(delta<0&&center+delta<lower)return (lower-center)/delta;
  return 1;
}
AMR_RECON_HD inline double amr_wall_pressure_value(int spatial_order,int reconstruction,
                                                    double probe_blend,const double* primitive,
                                                    const double* slopes,size_t n,size_t wall_id,
                                                    const AmrCell2D& cell,const WallFragment2D& wall,
                                                    const AmrWallPressureSample2D* samples) {
  const size_t owner=wall.owner_cell;
  const double center=primitive[3*n+owner];
  const double baseline=spatial_order==2
    ?amr_point_value(primitive,slopes,n,3,owner,cell,amr_wall_x(wall),amr_wall_y(wall)):center;
  if(reconstruction!=WALL_PRESSURE_NORMAL_PROBE||!(probe_blend>0)||!samples)return baseline;
  const auto sample=samples[wall_id];
  if(sample.first<0||sample.second<0||sample.first==sample.second||
     !(sample.first_distance>0)||!(sample.second_distance>sample.first_distance))return baseline;
  const double p1=primitive[3*n+static_cast<size_t>(sample.first)];
  const double p2=primitive[3*n+static_cast<size_t>(sample.second)];
  const double denominator=sample.second_distance-sample.first_distance;
  const double extrapolated=p1-sample.first_distance*(p2-p1)/denominator;
  double candidate=baseline+probe_blend*(extrapolated-baseline);
  if(!amr_vis_finite(candidate)||!(candidate>.1*center)||candidate>10*center)candidate=baseline;
  return candidate;
}
#undef AMR_RECON_HD
}
#endif
