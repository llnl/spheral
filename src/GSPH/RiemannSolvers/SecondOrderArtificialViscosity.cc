//---------------------------------Spheral++----------------------------------//
// SecondOrderArtificialViscosity
//   Frontiere, Raskin, Owen (2017) "CRKSPH:- A Conservative Reproducing Kernel 
//   Smoothed Particle Hydrodynamics Scheme," J. Comp. Phys.
//
// This is a reimplementation of the LimitedArtificialViscosity class as a
// derivative of RiemannSolverBase so it can be used with GSPH derived
// classes 
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#include "Hydro/HydroFieldNames.hh"
#include "GSPH/GSPHFieldNames.hh"

#include "GSPH/WaveSpeeds/WaveSpeedBase.hh"
#include "GSPH/Limiters/LimiterBase.hh"
#include "GSPH/RiemannSolvers/SecondOrderArtificialViscosity.hh"

#include <limits>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
SecondOrderArtificialViscosity<Dimension>::
SecondOrderArtificialViscosity(const Scalar Cl,
                               const Scalar Cq,
                               LimiterBase<Dimension>& slopeLimiter,
                               WaveSpeedBase<Dimension>& waveSpeed,
                               const bool linearReconstruction):
  RiemannSolverBase<Dimension>(slopeLimiter,
                               waveSpeed,
                               linearReconstruction),
  mCl(Cl),
  mCq(Cq){

}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SecondOrderArtificialViscosity<Dimension>::
~SecondOrderArtificialViscosity(){}


//------------------------------------------------------------------------------
// Interface State fluid hydro
//------------------------------------------------------------------------------
template<typename Dimension>
void
SecondOrderArtificialViscosity<Dimension>::
interfaceState(const typename Dimension::Vector& ri,
               const typename Dimension::Vector& rj,
               const typename Dimension::SymTensor& Hi,
               const typename Dimension::SymTensor& Hj,
               const typename Dimension::Scalar& rhoi,   
               const typename Dimension::Scalar& rhoj,
               const typename Dimension::Scalar& epsi,
               const typename Dimension::Scalar& epsj,
               const typename Dimension::Scalar& ci,   
               const typename Dimension::Scalar& cj,
               const typename Dimension::Scalar& Pi,    
               const typename Dimension::Scalar& Pj,
               const typename Dimension::Vector& vi,    
               const typename Dimension::Vector& vj,
               const typename Dimension::Vector& /*DrhoDxi*/,
               const typename Dimension::Vector& /*DrhoDxj*/,
               const typename Dimension::Vector& /*DepsDxi*/,
               const typename Dimension::Vector& /*DepsDxj*/,
               const typename Dimension::Vector& /*DpDxi*/,
               const typename Dimension::Vector& /*DpDxj*/,
               const typename Dimension::Tensor& DvDxi,
               const typename Dimension::Tensor& DvDxj,
                     typename Dimension::Scalar& Pstar,
                     typename Dimension::Vector& vstar,
                     typename Dimension::Scalar& /*rhostari*/,
                     typename Dimension::Scalar& /*rhostarj*/,
                     typename Dimension::Scalar& /*epsstari*/,
                     typename Dimension::Scalar& /*epsstarj*/) const{

  
  const auto tiny = std::numeric_limits<Scalar>::epsilon();

  const Vector rij = ri - rj;
  const Vector etaij = 0.5*(Hi+Hj)*rij;

  // default to nodal values
  Vector v1i = vi;
  Vector v1j = vj;

  // linear reconstruction
  if(this->linearReconstruction()){

    this->linearReconstruction(ri,rj, vi,vj,DvDxi,DvDxj, //inputs
                               v1i,v1j);                 //outputs
  
  }
  const Vector vij = v1i-v1j;
  const Scalar muij = std::max(0.0,-vij.dot(etaij)/(etaij.magnitude2() + tiny));
  const Scalar cij = 0.5*(ci+cj);
  const Scalar rhoij = 2*rhoi*rhoj/(rhoi+rhoj);
  Pstar = 0.5*(Pi+Pj) 
        + rhoij*muij*(this->Cl()*cij
                     +this->Cq()*muij);
  vstar = 0.5*(vi+vj);

}// Scalar interface class


template<typename Dimension>
void
SecondOrderArtificialViscosity<Dimension>::
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
