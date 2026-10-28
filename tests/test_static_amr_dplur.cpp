#include "static_amr2d_dplur.h"
#include "static_amr2d_gpu.h"
#include "static_amr2d_viscous.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace cfd;
namespace {
int failed(const std::string& message){std::cerr<<message<<'\n';return 1;}
double maximum(const std::vector<double>& values){double result=0;for(double value:values)result=std::max(result,std::fabs(value));return result;}
double difference(const std::vector<double>& a,const std::vector<double>& b){double result=0;for(size_t i=0;i<a.size();++i)result=std::max(result,std::fabs(a[i]-b[i]));return result;}
}

int main(){
  StaticAmrOptions2D options;options.nx=8;options.ny=6;options.levels=1;
  options.xmin=-1;options.xmax=1;options.ymin=-.75;options.ymax=.75;
  options.agglomerate_small=true;options.boxes.push_back({-.5,.5,-.4,.4,1});
  StaticAmrMesh2D mesh;std::string error;
  if(static_amr2d_build(options,{}, {}, {},&mesh,&error))return failed(error);
  Params2D params=params2d_default();params.amr_spatial_order=2;
  for(auto& boundary:mesh.boundaries){boundary.kind=AMR_BC_TRANSMISSIVE;boundary.reference[0]=1;boundary.reference[1]=.2;boundary.reference[2]=0;boundary.reference[3]=1/params.gamma;}
  if(static_amr2d_prepare_viscous(&mesh,params,&error))return failed(error);
  StaticAmrState2D uniform;const std::array<double,4> reference={1,.2,0,1/params.gamma};
  if(static_amr2d_initialize(mesh,params,reference,&uniform))return failed("uniform initialization failed");
  const auto original=uniform.u;StaticAmrDplurOptions2D controls;controls.iterations=5;controls.cfl=4;
  StaticAmrDplurReport2D report;
  if(static_amr2d_run_dplur(mesh,params,controls,&uniform,&report,&error))return failed(error);
  if(uniform.u!=original||report.accepted!=5||maximum(std::vector<double>(report.residual,report.residual+4))>2e-14)
    return failed("DP-LUR did not preserve uniform flow exactly");

  StaticAmrState2D full;if(static_amr2d_initialize(mesh,params,reference,&full))return 1;
  const size_t n=mesh.cells.size();const size_t center=n/2;
  full.u[3*n+center]*=1.02;
  std::vector<double> initial_rhs;if(static_amr2d_rhs(mesh,params,full.u,&initial_rhs))return failed("initial residual failed");
  const double initial=maximum(initial_rhs);controls.iterations=20;controls.cfl=.5;
  if(static_amr2d_run_dplur(mesh,params,controls,&full,&report,&error))return failed(error);
  std::vector<double> final_rhs;if(static_amr2d_rhs(mesh,params,full.u,&final_rhs))return failed("final residual failed");
  if(!(maximum(final_rhs)<initial))return failed("DP-LUR did not damp a pressure perturbation");

  StaticAmrState2D split;if(static_amr2d_initialize(mesh,params,reference,&split))return 1;
  split.u[3*n+center]*=1.02;controls.iterations=10;
  if(static_amr2d_run_dplur(mesh,params,controls,&split,&report,&error))return failed(error);
  if(static_amr2d_run_dplur(mesh,params,controls,&split,&report,&error))return failed(error);
  if(split.u!=full.u)return failed("DP-LUR 20-step and 10+10 restart states differ");
#ifdef CFD_HAVE_CUDA
  AmrGpuOptions gpu_options;gpu_options.dplur=true;gpu_options.dplur_sweeps=controls.sweeps;
  gpu_options.dplur_omega=controls.omega;gpu_options.dplur_block_coupling=controls.block_coupling;
  StaticAmrGpu* gpu=nullptr;StaticAmrState2D gpu_uniform;
  if(static_amr2d_initialize(mesh,params,reference,&gpu_uniform))return 1;
  if(static_amr_gpu_create(mesh,params,gpu_uniform,gpu_options,&gpu,&error))return failed(error);
  controls.iterations=5;controls.cfl=4;
  if(static_amr_gpu_run_dplur(gpu,controls,&gpu_uniform,&report,&error)){static_amr_gpu_destroy(gpu);return failed(error);}
  static_amr_gpu_destroy(gpu);gpu=nullptr;
  if(difference(gpu_uniform.u,original)>2e-14||report.accepted!=5)
    return failed("CUDA DP-LUR did not preserve uniform flow");

  StaticAmrState2D gpu_full;if(static_amr2d_initialize(mesh,params,reference,&gpu_full))return 1;
  gpu_full.u[3*n+center]*=1.02;controls.iterations=20;controls.cfl=.5;
  if(static_amr_gpu_create(mesh,params,gpu_full,gpu_options,&gpu,&error))return failed(error);
  if(static_amr_gpu_run_dplur(gpu,controls,&gpu_full,&report,&error)){static_amr_gpu_destroy(gpu);return failed(error);}
  static_amr_gpu_destroy(gpu);gpu=nullptr;
  if(difference(gpu_full.u,full.u)>2e-12)return failed("CPU/CUDA DP-LUR parity failed");

  StaticAmrState2D gpu_split;if(static_amr2d_initialize(mesh,params,reference,&gpu_split))return 1;
  gpu_split.u[3*n+center]*=1.02;controls.iterations=10;
  if(static_amr_gpu_create(mesh,params,gpu_split,gpu_options,&gpu,&error))return failed(error);
  if(static_amr_gpu_run_dplur(gpu,controls,&gpu_split,&report,&error)||
     static_amr_gpu_run_dplur(gpu,controls,&gpu_split,&report,&error)){static_amr_gpu_destroy(gpu);return failed(error);}
  static_amr_gpu_destroy(gpu);
  if(difference(gpu_split.u,gpu_full.u)>2e-14)
    return failed("CUDA DP-LUR 20-step and 10+10 restart states differ");
#endif
  return 0;
}
