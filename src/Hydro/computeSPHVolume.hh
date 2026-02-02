//---------------------------------Spheral++----------------------------------//
// Compute the volume from m/rho
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#ifndef __Spheral__computeSPHVolume__
#define __Spheral__computeSPHVolume__

#include <vector>

namespace Spheral {

// Forward declarations.
template<typename Dimension, typename DataType> class FieldList;


template<typename Dimension>
void
computeSPHVolume(const FieldList<Dimension, typename Dimension::Scalar>& mass,
                 const FieldList<Dimension, typename Dimension::Scalar>& massDensity,
                       FieldList<Dimension, typename Dimension::Scalar>& volume);

}

 
 #endif
