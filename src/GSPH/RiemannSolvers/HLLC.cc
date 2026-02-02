//---------------------------------Spheral++----------------------------------//
// HLLC -- approximate riemann solver
//    Toro E.F., Spruce M., Speares W., (1994) "Restoration of the Contact
//    Surface in the HLL-Riemann Solver," Shock Waves, 4:25-34
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#include "Hydro/HydroFieldNames.hh"
#include "GSPH/GSPHFieldNames.hh"

#include "GSPH/WaveSpeeds/WaveSpeedBase.hh"
#include "GSPH/Limiters/LimiterBase.hh"
#include "GSPH/RiemannSolvers/HLLC.hh"

#include <limits>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
HLLC<Dimension>::
HLLC(LimiterBase<Dimension>& slopeLimiter,
     WaveSpeedBase<Dimension>& waveSpeed,
     const bool linearReconstruction):
  RiemannSolverBase<Dimension>(slopeLimiter,
                               waveSpeed,
                               linearReconstruction){

}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
HLLC<Dimension>::
~HLLC(){}


//------------------------------------------------------------------------------
// Interface State fluid hydro
//------------------------------------------------------------------------------
template<typename Dimension>
void
HLLC<Dimension>::
interfaceState(const typename Dimension::Vector& ri,
               const typename Dimension::Vector& rj,
               const typename Dimension::SymTensor& Hi,
               const typename Dimension::SymTensor& Hj,
               const typename Dimension::Scalar& rhoi,   
               const typename Dimension::Scalar& rhoj, 
               const typename Dimension::Scalar& ci,   
               const typename Dimension::Scalar& cj,
               const typename Dimension::Scalar& Pi,    
               const typename Dimension::Scalar& Pj,
               const typename Dimension::Vector& vi,    
               const typename Dimension::Vector& vj,
               const typename Dimension::Vector& DrhoDxi,
               const typename Dimension::Vector& DrhoDxj,
               const typename Dimension::Vector& DpDxi,
               const typename Dimension::Vector& DpDxj,
               const typename Dimension::Tensor& DvDxi,
               const typename Dimension::Tensor& DvDxj,
                     typename Dimension::Scalar& Pstar,
                     typename Dimension::Vector& vstar,
                     typename Dimension::Scalar& rhostari,
                     typename Dimension::Scalar& rhostarj) const{

  Scalar Si, Sj;

  const auto tiny = std::numeric_limits<Scalar>::epsilon();

  const auto& waveSpeedObject = this->waveSpeed();

  const auto rij = ri - rj;
  const auto rhatij = rij.unitVector();

  vstar = 0.5*(vi+vj);
  Pstar = 0.5*(Pi+Pj);
  rhostari = rhoi;
  rhostarj = rhoj;

  if (ci > tiny or cj > tiny){


    // default to nodal values
    auto v1i = vi;
    auto v1j = vj;

    auto p1i = Pi;
    auto p1j = Pj;

    //auto rho1i = rhoi;
    //auto rho1j = rhoj;

    // linear reconstruction
    if(this->linearReconstruction()){

      // gradients along line of action
      //this->linearReconstruction(ri,rj, rhoi,rhoj,DrhoDxi,DrhoDxj, //inputs
      //                           rho1i,rho1j);                     //outputs
      this->linearReconstruction(ri,rj, Pi,Pj,DpDxi,DpDxj,         //inputs
                                 p1i,p1j);                         //outputs
      this->linearReconstruction(ri,rj, vi,vj,DvDxi,DvDxj,         //inputs
                                 v1i,v1j);                         //outputs
  
    }

    const auto ui = v1i.dot(rhatij);
    const auto uj = v1j.dot(rhatij);
    const auto wi = v1i - ui*rhatij;
    const auto wj = v1j - uj*rhatij;

    waveSpeedObject.waveSpeed(rhoi,rhoj,ci,cj,ui,uj,  //inputs
                              Si,Sj);                   //outputs

    const auto denom = safeInv(Si - Sj);

    const auto ustar = (Si*ui - Sj*uj - p1i + p1j )*denom;
    const auto wstar = (Si*wi - Sj*wj)*denom;
    vstar = ustar*rhatij + wstar;
    Pstar = Sj * (ustar-uj) + p1j;
    //rhostari = rho1i;// * (Si - ui)*safeInv(Si-ustar);
    //rhostarj = rho1j;// * (Sj - uj)*safeInv(Sj-ustar);

  }else{ // if ci & cj too small punt to normal av
    const auto uij = std::min((vi-vj).dot(rhatij),0.0);
    Pstar += 0.25 * (rhoi+rhoj) * (uij*uij);
  }
  
}// Scalar interface class


template<typename Dimension>
void
HLLC<Dimension>::
interfaceState(const Vector& /*ri*/,
               const Vector& /*rj*/,
               const SymTensor& /*Hi*/,
               const SymTensor& /*Hj*/,
               const Scalar& /*rhoi*/,   
               const Scalar& /*rhoj*/, 
               const Scalar& /*ci*/,   
               const Scalar& /*cj*/,
               const Scalar& /*Pi*/,    
               const Scalar& /*Pj*/,
               const Vector& /*vi*/,    
               const Vector& /*vj*/,
               const SymTensor& /*Si*/,    
               const SymTensor& /*Sj*/,
               const Tensor& /*Di*/,    
               const Tensor& /*Dj*/,
                     Vector& /*Tstar*/,
                     Vector& /*vstar*/) const{



}

} // spheral namespace