//---------------------------------Spheral++----------------------------------//
// EinfeldtWaveSpeed 
//   Einfeldt B., "On Godunov-Type Methods for Gas Dynamics," Siam. J. Numer.
//   Anal., 25(2):294-318, 1988
//   
//   Nico Fleischmann, Stefan Adami, Nikolaus A. Adams, On an inconsistency of the 
//   arithmetic-average signal speed estimate for HLL-type Riemann solvers, Journal of 
//   Computational Physics: X, Volume 8, 2020.
//
//   Fleischmann reccomended modifying the Einfeldt wave speed with the min/max
//   treatment like Davis. Also compared the accuracy of the two wave speeds.
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#include "Utilities/safeInv.hh"
#include "EinfeldtWaveSpeed.hh"
#include <math.h>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
EinfeldtWaveSpeed<Dimension>::
EinfeldtWaveSpeed(){}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
EinfeldtWaveSpeed<Dimension>::
~EinfeldtWaveSpeed(){}

//------------------------------------------------------------------------------
// wave speed
//------------------------------------------------------------------------------
template<typename Dimension>
void
EinfeldtWaveSpeed<Dimension>::
waveSpeed(const typename Dimension::Scalar rhoi, 
          const typename Dimension::Scalar rhoj, 
          const typename Dimension::Scalar ci,   
          const typename Dimension::Scalar cj, 
          const typename Dimension::Scalar ui,   
          const typename Dimension::Scalar uj,
                typename Dimension::Scalar& Si,  
                typename Dimension::Scalar& Sj) const {

  const auto tiny = std::numeric_limits<Scalar>::epsilon();
  
  const auto sRhoi = sqrt(std::max(rhoi,tiny));
  const auto sRhoj = sqrt(std::max(rhoj,tiny));
  const auto denom = 1.0/std::max(sRhoi+sRhoj,tiny);
  
  const auto eta = 0.5*sRhoi*sRhoj*denom*denom;
  
  const auto d2 = (sRhoi * ci*ci + sRhoj * cj*cj)*denom + eta * (ui-uj)*(ui-uj);

  const auto d = sqrt(d2);
  const auto utilde = (sRhoi * ui + sRhoj * uj)*denom;

  
  Si = std::max(utilde + d, ui + ci);
  Sj = std::min(utilde - d, uj - cj);
}


}