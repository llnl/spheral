//---------------------------------Spheral++----------------------------------//
// AcousticWaveSpeed -- simple acoustic wave speed
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#include "AcousticWaveSpeed.hh"
#include <stdio.h>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
AcousticWaveSpeed<Dimension>::
AcousticWaveSpeed(){}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
AcousticWaveSpeed<Dimension>::
~AcousticWaveSpeed(){}

//------------------------------------------------------------------------------
// wave speed
//------------------------------------------------------------------------------
template<typename Dimension>
void
AcousticWaveSpeed<Dimension>::
waveSpeed(const typename Dimension::Scalar /*rhoi*/, 
          const typename Dimension::Scalar /*rhoj*/, 
          const typename Dimension::Scalar ci,   
          const typename Dimension::Scalar cj, 
          const typename Dimension::Scalar ui,   
          const typename Dimension::Scalar uj,
                typename Dimension::Scalar& Si,  
                typename Dimension::Scalar& Sj) const {
  Si = ui+ci;
  Sj = uj-cj;
}


}