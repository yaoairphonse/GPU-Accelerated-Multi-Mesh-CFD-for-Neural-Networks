#include "static_amr2d_diagnostics.h"

#include "static_amr2d_sa.h"
#include "static_amr2d_viscous.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>

namespace cfd {
namespace {
struct Budget {double l1[5]={},signed_value[5]={};};
bool failed(std::string* error,const std::string& message){if(error)*error=message;return false;}
void add(Budget& budget,const double* value,int nv){for(int v=0;v<nv;++v){budget.l1[v]+=std::fabs(value[v]);budget.signed_value[v]+=value[v];}}
double drag(double fx,double fy,double angle){return fx*std::cos(angle)+fy*std::sin(angle);}
double lift(double fx,double fy,double angle){return-fx*std::sin(angle)+fy*std::cos(angle);}
}

int static_amr2d_write_budgets(const StaticAmrMesh2D& mesh,const Params2D& p,
                               const StaticAmrState2D& state,
                               const StaticAmrAuditOptions2D& options,
                               const std::string& prefix,std::string* error){
  const size_t n=mesh.cells.size(),nf=mesh.faces.size(),nw=mesh.walls.size();
  const bool use_sa=mesh.transport.sa_model;const int nv=use_sa?5:4;
  if(prefix.empty()||n==0||state.u.size()!=static_cast<size_t>(nv)*n||
     !(options.dynamic_pressure>0)||!(options.reference_length>0)||
     options.bins<1||options.bins>10000){failed(error,"invalid 2D audit controls");return 1;}
  std::vector<double> residual,face_flux;
  if(static_amr2d_rhs(mesh,p,state.u,&residual,&face_flux)){failed(error,"2D audit residual failed");return 1;}
  std::vector<double> primitive(4*n);
  for(size_t k=0;k<n;++k){const double rho=state.u[k],u=state.u[n+k]/rho,v=state.u[2*n+k]/rho;
    primitive[k]=rho;primitive[n+k]=u;primitive[2*n+k]=v;
    primitive[3*n+k]=(p.gamma-1)*(state.u[3*n+k]-.5*rho*(u*u+v*v));}
  AmrViscousFields2D viscous;SaFields2D sa;
  if(mesh.transport.model==AMR_MODEL_LAMINAR&&static_amr2d_viscous_evaluate(mesh,p,primitive,&viscous,&state.u)){
    failed(error,"2D audit viscous evaluation failed");return 1;}
  if(use_sa&&static_amr2d_sa_evaluate(mesh,p,state.u,primitive,viscous,face_flux,&sa)){
    failed(error,"2D audit SA evaluation failed");return 1;}
  Budget same_level,two_to_one,farfield,wall,source;
  for(size_t q=0;q<nf;++q){const auto& face=mesh.faces[q];double integrated[5]={};
    for(int v=0;v<4;++v)integrated[v]=face_flux[v*nf+q];
    if(use_sa)integrated[4]=sa.face_flux[q];
    Budget* family=face.boundary>=0?&farfield:
      (mesh.cells[face.left].level==mesh.cells[face.right].level?&same_level:&two_to_one);
    if(face.left>=0){double contribution[5];for(int v=0;v<nv;++v)contribution[v]=-integrated[v];add(*family,contribution,nv);}
    if(face.right>=0){double contribution[5];for(int v=0;v<nv;++v)contribution[v]=integrated[v];add(*family,contribution,nv);}
  }
  std::vector<double> wall_pressure;if(static_amr2d_wall_pressure(mesh,p,state.u,&wall_pressure)){
    failed(error,"2D audit wall pressure failed");return 1;}
  struct Load {long fragments=0;double pfx=0,pfy=0,vfx=0,vfy=0,pm=0,vm=0;};
  double xmin=std::numeric_limits<double>::infinity(),xmax=-xmin;
  for(const auto& fragment:mesh.walls){const double x=.5*(fragment.x0+fragment.x1);xmin=std::min(xmin,x);xmax=std::max(xmax,x);}
  std::vector<Load> bins(options.bins);
  for(size_t q=0;q<nw;++q){const auto& fragment=mesh.walls[q];const double x=.5*(fragment.x0+fragment.x1),y=.5*(fragment.y0+fragment.y1);
    const double pfx=wall_pressure[q]*fragment.nx*fragment.length,pfy=wall_pressure[q]*fragment.ny*fragment.length;
    const double vfx=mesh.transport.model==AMR_MODEL_LAMINAR?-viscous.wall_flux[nw+q]:0;
    const double vfy=mesh.transport.model==AMR_MODEL_LAMINAR?-viscous.wall_flux[2*nw+q]:0;
    double contribution[5]={};
    // Conserved residual ordering is rho, rho*u, rho*v, rho*E, rho*nutilde.
    contribution[0]=0;contribution[1]=-pfx-vfx;contribution[2]=-pfy-vfy;contribution[3]=0;
    if(use_sa)contribution[4]=sa.wall_flux[q];add(wall,contribution,nv);
    const double fraction=xmax>xmin?(x-xmin)/(xmax-xmin):0;
    const int index=std::min(options.bins-1,std::max(0,static_cast<int>(fraction*options.bins)));
    auto& bin=bins[index];++bin.fragments;bin.pfx+=pfx;bin.pfy+=pfy;bin.vfx+=vfx;bin.vfy+=vfy;
    bin.pm+=(x-options.moment_x)*pfy-(y-options.moment_y)*pfx;
    bin.vm+=(x-options.moment_x)*vfy-(y-options.moment_y)*vfx;
  }
  if(use_sa)for(size_t k=0;k<n;++k){const double integrated=mesh.cells[k].volume*sa.source[k];source.l1[4]+=std::fabs(integrated);source.signed_value[4]+=integrated;}
  std::ofstream face_file(prefix+".face_budget.csv");
  face_file<<"family,l1_rho,l1_mx,l1_my,l1_E,l1_SA,signed_integral_rho,signed_integral_mx,signed_integral_my,signed_integral_E,signed_integral_SA\n";
  const auto row=[&](const char* name,const Budget& budget){face_file<<name;for(double value:budget.l1)face_file<<','<<std::setprecision(17)<<value;for(double value:budget.signed_value)face_file<<','<<value;face_file<<'\n';};
  row("interface_same_level",same_level);row("interface_2to1",two_to_one);row("farfield",farfield);row("wall",wall);row("SA_source",source);face_file.close();
  if(!face_file){failed(error,"2D face budget write failed");return 1;}
  const double qref=options.dynamic_pressure,cmref=qref*options.reference_length;
  std::ofstream load_file(prefix+".wall_load_budget.csv");
  load_file<<"axis,bin,coordinate_min,coordinate_max,fragments,CD_pressure,CD_viscous,CL_pressure,CL_viscous,CM_z_pressure,CM_z_viscous\n";
  const double width=xmax>xmin?(xmax-xmin)/options.bins:0;
  for(int index=0;index<options.bins;++index){const auto& bin=bins[index];load_file<<"x,"<<index<<','<<std::setprecision(17)<<xmin+index*width<<','<<xmin+(index+1)*width<<','<<bin.fragments<<','
      <<drag(bin.pfx,bin.pfy,options.angle_radians)/qref<<','<<drag(bin.vfx,bin.vfy,options.angle_radians)/qref<<','
      <<lift(bin.pfx,bin.pfy,options.angle_radians)/qref<<','<<lift(bin.vfx,bin.vfy,options.angle_radians)/qref<<','
      <<bin.pm/cmref<<','<<bin.vm/cmref<<'\n';}
  load_file.close();if(!load_file){failed(error,"2D wall-load budget write failed");return 1;}
  return 0;
}
} // namespace cfd
