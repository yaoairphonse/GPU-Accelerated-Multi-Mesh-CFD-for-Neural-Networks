#include "static_amr2d.h"
#include "static_amr2d_gpu.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace cfd;

namespace {
double maximum(const std::vector<double>& values) {
  double result=0;
  for(double value:values)result=std::max(result,std::fabs(value));
  return result;
}

int fail(const std::string& message) {
  std::cerr<<message<<'\n';
  return 1;
}
}

int main() {
  StaticAmrOptions2D options;
  options.nx=8;options.ny=8;options.levels=2;options.buffer=4;
  options.xmin=-2;options.xmax=2;options.ymin=-2;options.ymax=2;
  options.max_cells=100000;
  const std::vector<double> x={-.61,.04,.57,-.08};
  const std::vector<double> y={-.07,-.53,.11,.59};
  const std::vector<int> patches(x.size(),1);
  StaticAmrMesh2D mesh;std::string error;
  if(static_amr2d_build(options,x,y,patches,&mesh,&error))return fail(error);

  Params2D params=params2d_default();params.cfl=.1;params.amr_spatial_order=2;
  params.wall_pressure_reconstruction=WALL_PRESSURE_NORMAL_PROBE;
  params.wall_pressure_probe_blend=1;
  if(static_amr2d_prepare_viscous(&mesh,params,&error))return fail(error);
  if(mesh.wall_pressure_samples.size()!=mesh.walls.size())return fail("wall-pressure donor table is incomplete");
  for(const auto& sample:mesh.wall_pressure_samples)if(sample.first<0||sample.second<0||sample.first==sample.second||
    !(sample.first_distance>0)||!(sample.second_distance>sample.first_distance))return fail("invalid wall-pressure donor geometry");

  const std::array<double,4> rest={1,0,0,1.0/params.gamma};
  StaticAmrState2D state;if(static_amr2d_initialize(mesh,params,rest,&state))return fail("initialization failed");
  std::vector<double> rhs,pressure;
  if(static_amr2d_rhs(mesh,params,state.u,&rhs)||maximum(rhs)>2e-11)return fail("normal-probe uniform-rest residual is not zero");
  if(static_amr2d_wall_pressure(mesh,params,state.u,&pressure))return fail("wall-pressure evaluation failed");
  for(double value:pressure)if(std::fabs(value-rest[3])>2e-14)return fail("normal-probe constant pressure is not exact");

  StaticAmrState2D nonlinear=state;const size_t cell_count=mesh.cells.size();
  for(size_t k=0;k<cell_count;++k){const double value=.8+.07*mesh.cells[k].cx*mesh.cells[k].cx+.03*mesh.cells[k].cy;
    nonlinear.u[3*cell_count+k]=value/(params.gamma-1);}
  std::vector<double> probe_pressure;if(static_amr2d_wall_pressure(mesh,params,nonlinear.u,&probe_pressure))return fail("nonlinear probe pressure failed");
  for(double value:probe_pressure)if(!(value>0)||!std::isfinite(value))return fail("nonlinear probe pressure is not positive and finite");

  Params2D centroid=params2d_default();centroid.cfl=.1;centroid.amr_spatial_order=2;
  if(static_amr2d_prepare_viscous(&mesh,centroid,&error))return fail(error);
  std::vector<double> centroid_pressure;if(static_amr2d_wall_pressure(mesh,centroid,nonlinear.u,&centroid_pressure))return fail("nonlinear centroid pressure failed");
  double pressure_difference=0;for(size_t q=0;q<probe_pressure.size();++q)pressure_difference=std::max(pressure_difference,std::fabs(probe_pressure[q]-centroid_pressure[q]));
  if(!(pressure_difference>1e-8))return fail("normal-probe option did not change a nonlinear wall-pressure field");
  std::vector<double> baseline;if(static_amr2d_rhs(mesh,centroid,state.u,&baseline))return fail("centroid residual failed");
  Params2D zero_blend=centroid;zero_blend.wall_pressure_reconstruction=WALL_PRESSURE_NORMAL_PROBE;
  zero_blend.wall_pressure_probe_blend=0;
  if(static_amr2d_prepare_viscous(&mesh,zero_blend,&error))return fail(error);
  std::vector<double> disabled;if(static_amr2d_rhs(mesh,zero_blend,state.u,&disabled))return fail("zero-blend residual failed");
  if(baseline!=disabled)return fail("disabled normal-probe mode changed the baseline residual");

  char reason[256]={};
  if(solver2d_cuda_available(reason,sizeof(reason))){
    zero_blend.wall_pressure_probe_blend=1;
    if(static_amr2d_prepare_viscous(&mesh,zero_blend,&error))return fail(error);
    StaticAmrGpu* gpu=nullptr;AmrGpuOptions gpu_options;gpu_options.batch_trials=1;gpu_options.graph=false;
    if(static_amr_gpu_create(mesh,zero_blend,state,gpu_options,&gpu,&error))return fail(error);
    const int status=static_amr_gpu_advance(gpu,1e-4,100000,&state,&error);
    if(!status&&static_amr_gpu_download(gpu,&state,&error))return fail(error);
    static_amr_gpu_destroy(gpu);
    if(status)return fail(error);
    const size_t n=mesh.cells.size();
    for(size_t k=0;k<n;++k)if(std::fabs(state.u[k]-1)>2e-11||std::fabs(state.u[n+k])>2e-11||
      std::fabs(state.u[2*n+k])>2e-11||std::fabs(state.u[3*n+k]-rest[3]/(zero_blend.gamma-1))>2e-11)
      return fail("CUDA normal-probe uniform-rest state changed");
  }
  return 0;
}
