//---------------------------------Spheral++----------------------------------//
// RiemannSolverBase
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Boundary/Boundary.hh"

#include "Hydro/HydroFieldNames.hh"
#include "GSPH/GSPHFieldNames.hh"

#include "GSPH/WaveSpeeds/WaveSpeedBase.hh"
#include "GSPH/Limiters/LimiterBase.hh"
#include "GSPH/RiemannSolvers/RiemannSolverBase.hh"


#ifdef _OPENMP
#include "omp.h"
#endif

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
RiemannSolverBase<Dimension>::
RiemannSolverBase(LimiterBase<Dimension>& slopeLimiter,
                  WaveSpeedBase<Dimension>& waveSpeed,
                  const bool linearReconstruction):
  mSlopeLimiter(slopeLimiter),
  mWaveSpeed(waveSpeed),
  mLinearReconstruction(linearReconstruction){
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
RiemannSolverBase<Dimension>::
~RiemannSolverBase(){}

//------------------------------------------------------------------------------
// non-op
//------------------------------------------------------------------------------
template<typename Dimension>
void
RiemannSolverBase<Dimension>::
initialize(const DataBase<Dimension>& /*dataBase*/,
           const State<Dimension>& /*state*/,
           const StateDerivatives<Dimension>& /*derivs*/,
           typename RiemannSolverBase<Dimension>::ConstBoundaryIterator /*boundaryBegin*/,
           typename RiemannSolverBase<Dimension>::ConstBoundaryIterator /*boundaryEnd*/,
           const typename Dimension::Scalar /*time*/,
           const typename Dimension::Scalar /*dt*/,
           const TableKernel<Dimension>& /*W*/){

} // initialize method


//------------------------------------------------------------------------------
// reconstruct from limited gradient
//------------------------------------------------------------------------------
template<typename Dimension>
void
RiemannSolverBase<Dimension>::
linearReconstruction(const typename Dimension::Vector& ri,
                     const typename Dimension::Vector& rj,
                     const typename Dimension::Scalar& yi,
                     const typename Dimension::Scalar& yj,
                     const typename Dimension::Vector& DyDxi,
                     const typename Dimension::Vector& DyDxj,
                           typename Dimension::Scalar& ytildei,
                           typename Dimension::Scalar& ytildej) const {
  
  const auto tiny = std::numeric_limits<Scalar>::epsilon();

  const auto rij = (ri-rj);

  // relavant deltas in field value
  const auto Dy0 = (yi-yj);
  const auto Dyi = 0.5*DyDxi.dot(rij);
  const auto Dyj = 0.5*DyDxj.dot(rij);

  // ratios of SPH derivs to ij particle difference
  const auto denom = 2.0 / (sgn(Dy0) * std::max(tiny,abs(Dy0)));
  const auto ratioi = Dyi * denom;
  const auto ratioj = Dyj * denom;

  // limiter function 
  const auto phii = this->mSlopeLimiter.slopeLimiter(ratioi);
  const auto phij = this->mSlopeLimiter.slopeLimiter(ratioj);
  const auto phi = std::min(phii,phij);

  // linear constructed inteface values
  ytildei = yi - phi * Dyi;
  ytildej = yj + phi * Dyj;
}

//------------------------------------------------------------------------------
// reconstruct from limited gradient
//------------------------------------------------------------------------------
template<typename Dimension>
void
RiemannSolverBase<Dimension>::
linearReconstruction(const typename Dimension::Vector& ri,
                     const typename Dimension::Vector& rj,
                     const typename Dimension::Vector& yi,
                     const typename Dimension::Vector& yj,
                     const typename Dimension::Tensor& DyDxi,
                     const typename Dimension::Tensor& DyDxj,
                           typename Dimension::Vector& ytildei,
                           typename Dimension::Vector& ytildej) const {
  
  const auto tiny = std::numeric_limits<Scalar>::epsilon();

  const auto rij = (ri-rj);
  const auto rhatij = rij.unitVector();

  // relavant deltas in field value
  const auto Dy0 = (yi-yj);
  const auto Dyi = 0.5*DyDxi.dot(rij);   
  const auto Dyj = 0.5*DyDxj.dot(rij);  

  // scalar delta along line of action s
  const auto Dy0s = Dy0.dot(rhatij);
  const auto Dyis = Dyi.dot(rhatij);
  const auto Dyjs = Dyj.dot(rhatij);

  // ratios of SPH derivs to ij particle difference
  const auto denom = 2.0 / (sgn(Dy0s) * std::max(tiny,abs(Dy0s)));
  const auto ratioi = Dyis * denom;
  const auto ratioj = Dyjs * denom;

  // limiter function 
  const auto phii = this->mSlopeLimiter.slopeLimiter(ratioi);
  const auto phij = this->mSlopeLimiter.slopeLimiter(ratioj);
  const auto phi = std::min(phii,phij);

  ytildei = yi - phi * Dyi;
  ytildej = yj + phi * Dyj;
}

//------------------------------------------------------------------------------
// default to non-op
//------------------------------------------------------------------------------


template<typename Dimension>
void
RiemannSolverBase<Dimension>::
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
               const Vector& /*DrhoDxi*/,    
               const Vector& /*DrhoDxj*/,
               const Vector& /*DpDxi*/,    
               const Vector& /*DpDxj*/,
               const Tensor& /*DvDxi*/,    
               const Tensor& /*DvDxj*/,
                     Scalar& /*Pstar*/,
                     Vector& /*vstar*/,
                     Scalar& /*rhostari*/,
                     Scalar& /*rhostarj*/) const{

}

//------------------------------------------------------------------------------
// default to non-op
//------------------------------------------------------------------------------
template<typename Dimension>
void
RiemannSolverBase<Dimension>::
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