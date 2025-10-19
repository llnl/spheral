//---------------------------------Spheral++----------------------------------//
// DavisWaveSpeed 
//   Davis S., (1988) "Simplified Second-Order Godunov-Type Methods" Siam. Sci. 
//   Stat. Comput., 9(3):445-473
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#include <algorithm>

#include "DavisWaveSpeed.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
DavisWaveSpeed<Dimension>::
DavisWaveSpeed(){}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
DavisWaveSpeed<Dimension>::
~DavisWaveSpeed(){}

//------------------------------------------------------------------------------
// wave speed
//------------------------------------------------------------------------------
template<typename Dimension>
void
DavisWaveSpeed<Dimension>::
waveSpeed(const typename Dimension::Scalar rhoi, 
          const typename Dimension::Scalar rhoj, 
          const typename Dimension::Scalar ci,   
          const typename Dimension::Scalar cj, 
          const typename Dimension::Scalar ui,   
          const typename Dimension::Scalar uj,
                typename Dimension::Scalar& Si,  
                typename Dimension::Scalar& Sj) const {
  Si = std::max(uj+cj,ui+ci);
  Sj = std::min(ui-ci,uj-cj);
}


}