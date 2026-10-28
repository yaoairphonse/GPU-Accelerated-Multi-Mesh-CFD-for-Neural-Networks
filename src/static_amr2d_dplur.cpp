// Matrix-free block DP-LUR reference path for steady 2-D Cartesian AMR.
#include "static_amr2d_dplur.h"

#include "low_mach.h"
#include "static_amr2d_sa.h"
#include "static_amr2d_viscous.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>

namespace cfd {
namespace {

bool fail(std::string* error,const std::string& message) {
  if(error)*error=message;
  return false;
}

bool primitive(const std::vector<double>& u,size_t n,double gamma,
               bool sa,std::vector<double>* w) {
  if(u.size()!=(sa?5:4)*n)return false;
  w->assign(4*n,0);
  for(size_t k=0;k<n;++k) {
    const double rho=u[k],mx=u[n+k],my=u[2*n+k],energy=u[3*n+k];
    if(!(rho>0)||!std::isfinite(rho)||!std::isfinite(mx)||
       !std::isfinite(my)||!std::isfinite(energy))return false;
    const double vx=mx/rho,vy=my/rho;
    const double pressure=(gamma-1)*(energy-.5*(mx*vx+my*vy));
    if(!(pressure>0)||!std::isfinite(pressure))return false;
    (*w)[k]=rho;(*w)[n+k]=vx;(*w)[2*n+k]=vy;(*w)[3*n+k]=pressure;
    if(sa&&!std::isfinite(u[4*n+k]))return false;
  }
  return true;
}

void euler_normal_jacobian(const double* w,const double* du,
                           double nx,double ny,double gamma,double* action) {
  const double rho=w[0],u=w[1],v=w[2],p=w[3],drho=du[0];
  const double un=u*nx+v*ny,speed2=u*u+v*v;
  const double energy=p/(gamma-1)+.5*rho*speed2;
  const double velocity_x=(du[1]-u*drho)/rho;
  const double velocity_y=(du[2]-v*drho)/rho;
  const double velocity_normal=velocity_x*nx+velocity_y*ny;
  const double pressure=(gamma-1)*(du[3]-u*du[1]-v*du[2]+.5*speed2*drho);
  action[0]=du[1]*nx+du[2]*ny;
  action[1]=du[1]*un+rho*u*velocity_normal+pressure*nx;
  action[2]=du[2]*un+rho*v*velocity_normal+pressure*ny;
  action[3]=(du[3]+pressure)*un+(energy+p)*velocity_normal;
}

double face_speed(const Params2D& p,const std::vector<double>& w,size_t n,
                  int cell,int axis) {
  const double velocity=w[(axis+1)*n+cell];
  const double sound=std::sqrt(p.gamma*w[3*n+cell]/w[cell]);
  if(!p.low_mach_preconditioning)return std::fabs(velocity)+sound;
  const double primitive[4]={w[cell],w[n+cell],w[2*n+cell],w[3*n+cell]};
  const double other[4]={primitive[0],primitive[1],primitive[2],primitive[3]};
  return lm_face_speed(primitive,other,axis+1,3,p.gamma,
                       p.low_mach_reference,p.low_mach_floor);
}

void residual_norm(const StaticAmrMesh2D& mesh,const std::vector<double>& r,
                   bool sa,double reference_sa,double* output) {
  const size_t n=mesh.cells.size();const int nv=sa?5:4;
  double volume=0;for(int v=0;v<5;++v)output[v]=0;
  for(size_t k=0;k<n;++k){const double cell_volume=mesh.cells[k].volume;
    volume+=cell_volume;for(int v=0;v<nv;++v)output[v]+=cell_volume*std::fabs(r[v*n+k]);}
  for(int v=0;v<nv;++v)output[v]/=volume;
  if(sa&&reference_sa>0)output[4]/=reference_sa;
}

} // namespace

int static_amr2d_run_dplur(const StaticAmrMesh2D& mesh,const Params2D& p,
                           const StaticAmrDplurOptions2D& requested,
                           StaticAmrState2D* state,StaticAmrDplurReport2D* report,
                           std::string* error) {
  const size_t n=mesh.cells.size();const bool use_sa=mesh.transport.sa_model;
  const int nv=use_sa?5:4;
  if(report)*report={};
  if(!state||n==0||state->u.size()!=static_cast<size_t>(nv)*n||
     requested.iterations<0||requested.sweeps<1||requested.sweeps>50||
     !(requested.cfl>0)||!std::isfinite(requested.cfl)||
     !(requested.omega>0&&requested.omega<=1)||
     !std::isfinite(requested.budget_seconds)||requested.budget_seconds<0||
     !(p.gamma>1)||!std::isfinite(p.gamma)||
     (p.flux!=FLUX_HLLC&&p.flux!=FLUX_AUSMPUP)||
     mesh.small_cells!=0||!mesh.options.agglomerate_small) {
    fail(error,"DP-LUR requires valid controls and an agglomerated AMR mesh without SRD cells");
    return 1;
  }
  const auto start=std::chrono::steady_clock::now();
  const auto elapsed=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
  std::vector<double> residual,face_flux,w,mean_rate(n),convective_rate(n),
    viscous_rate(n),sa_rate(n),dt(n),diagonal(static_cast<size_t>(nv)*n),
    update(static_cast<size_t>(nv)*n),next(update.size()),sum(update.size()),candidate(update.size());
  double cfl=requested.cfl;long accepted=0,rejected=0;

  while(accepted<requested.iterations) {
    if(requested.budget_seconds>0&&elapsed()>=requested.budget_seconds) {
      state->rejected+=rejected;
      if(report){report->accepted=accepted;report->rejected=rejected;report->seconds=elapsed();report->final_cfl=cfl;}
      return 2;
    }
    if(static_amr2d_rhs(mesh,p,state->u,&residual,&face_flux)||
       !primitive(state->u,n,p.gamma,use_sa,&w)) {
      fail(error,"DP-LUR base residual contains a nonphysical state");return 1;
    }
    if(p.low_mach_preconditioning)for(size_t k=0;k<n;++k){
      double rr[5]={},ww[5]={w[k],w[n+k],w[2*n+k],w[3*n+k],
                             use_sa?state->u[4*n+k]/w[k]:0};
      for(int v=0;v<nv;++v)rr[v]=residual[v*n+k];
      lm_apply(rr,ww,p.gamma,p.low_mach_reference,p.low_mach_floor,2,use_sa);
      for(int v=0;v<nv;++v)residual[v*n+k]=rr[v];
    }

    std::fill(convective_rate.begin(),convective_rate.end(),0);
    for(const auto& face:mesh.faces){const int inside=face.left>=0?face.left:face.right;
      double speed=face_speed(p,w,n,inside,face.axis);
      if(face.left>=0&&face.right>=0)speed=std::max(speed,face_speed(p,w,n,face.right,face.axis));
      const double integrated=speed*face.length;
      if(face.left>=0)convective_rate[face.left]+=integrated/mesh.cells[face.left].volume;
      if(face.right>=0)convective_rate[face.right]+=integrated/mesh.cells[face.right].volume;
    }
    for(const auto& wall:mesh.walls){const size_t k=wall.owner_cell;
      const double sound=std::sqrt(p.gamma*w[3*n+k]/w[k]);
      convective_rate[k]+=sound*wall.length/mesh.cells[k].volume;
    }
    std::fill(viscous_rate.begin(),viscous_rate.end(),0);
    std::fill(sa_rate.begin(),sa_rate.end(),0);
    if(mesh.transport.model==AMR_MODEL_LAMINAR){
      AmrViscousFields2D viscous;
      if(static_amr2d_viscous_evaluate(mesh,p,w,&viscous,&state->u)){
        fail(error,"DP-LUR viscous-rate evaluation failed");return 1;}
      for(size_t k=0;k<n;++k)if(std::isfinite(viscous.cell_dt[k])&&viscous.cell_dt[k]>0)
        viscous_rate[k]=mesh.transport.viscous_cfl/viscous.cell_dt[k];
      if(use_sa){SaFields2D sa;
        if(static_amr2d_sa_evaluate(mesh,p,state->u,w,viscous,face_flux,&sa)){
          fail(error,"DP-LUR SA-rate evaluation failed");return 1;}
        for(size_t k=0;k<n;++k)if(std::isfinite(sa.cell_dt[k])&&sa.cell_dt[k]>0)
          sa_rate[k]=mesh.transport.sa_cfl/sa.cell_dt[k];
      }
    }
    for(size_t k=0;k<n;++k){
      mean_rate[k]=convective_rate[k]+2*viscous_rate[k];
      const double rate=std::max(mean_rate[k],use_sa?sa_rate[k]:0.);
      if(!(rate>0)||!std::isfinite(rate)){fail(error,"DP-LUR local spectral rate is invalid");return 1;}
      dt[k]=cfl/rate;
      diagonal[k]=1+dt[k]*convective_rate[k];
      for(int v=1;v<4;++v)diagonal[v*n+k]=1+dt[k]*mean_rate[k];
      if(use_sa)diagonal[4*n+k]=1+dt[k]*std::max(sa_rate[k],convective_rate[k]);
      for(int v=0;v<nv;++v)update[v*n+k]=dt[k]*residual[v*n+k]/diagonal[v*n+k];
    }

    for(int sweep=1;sweep<requested.sweeps;++sweep){
      std::fill(sum.begin(),sum.end(),0);
      for(const auto& face:mesh.faces)if(face.left>=0&&face.right>=0){
        const int cells[2]={face.left,face.right};
        const double speed=std::max(face_speed(p,w,n,face.left,face.axis),
                                    face_speed(p,w,n,face.right,face.axis));
        for(int side=0;side<2;++side){const int cell=cells[side],neighbor=cells[side^1];
          const double inverse_volume=1/mesh.cells[cell].volume;
          const double convective=.5*speed*face.length*inverse_volume;
          double neighbor_update[4],action[4];
          for(int v=0;v<4;++v)neighbor_update[v]=update[v*n+neighbor];
          if(requested.block_coupling){
            const double normal=(side==0?1.:-1.);
            const double primitive_neighbor[4]={w[neighbor],w[n+neighbor],w[2*n+neighbor],w[3*n+neighbor]};
            euler_normal_jacobian(primitive_neighbor,neighbor_update,
                                  face.axis==0?normal:0,face.axis==1?normal:0,
                                  p.gamma,action);
            const double scale=.5*face.length*inverse_volume;
            for(int v=0;v<4;++v)sum[v*n+cell]+=convective*neighbor_update[v]-scale*action[v];
          }else for(int v=0;v<4;++v)sum[v*n+cell]+=convective*neighbor_update[v];
          if(use_sa)sum[4*n+cell]+=convective*update[4*n+neighbor];
        }
      }
      for(size_t k=0;k<n;++k)for(int v=0;v<nv;++v){
        const size_t q=v*n+k;
        const double value=dt[k]*(residual[q]+sum[q])/diagonal[q];
        next[q]=(1-requested.omega)*update[q]+requested.omega*value;
      }
      update.swap(next);
    }
    for(size_t q=0;q<candidate.size();++q)candidate[q]=state->u[q]+update[q];
    if(!primitive(candidate,n,p.gamma,use_sa,&w)){
      ++rejected;cfl*=.5;
      if(!(cfl>=1e-6)||rejected>20){fail(error,"DP-LUR candidates remained nonphysical after CFL reductions");return 1;}
      continue;
    }
    state->u.swap(candidate);++state->steps;++accepted;
    state->time=static_cast<double>(state->steps); // explicit pseudo-iteration coordinate
  }

  if(static_amr2d_rhs(mesh,p,state->u,&residual)){
    fail(error,"DP-LUR final residual evaluation failed");return 1;}
  state->rejected+=rejected;
  if(report){report->accepted=accepted;report->rejected=rejected;
    report->seconds=elapsed();report->final_cfl=cfl;
    residual_norm(mesh,residual,use_sa,mesh.transport.sa_nu_inf,report->residual);}
  return 0;
}

} // namespace cfd
