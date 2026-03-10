//---------------------------------Spheral++----------------------------------//
// planarReflectingOperator
// Generate the reflection operator for the given plane.
//
// Created by JMO, Wed Feb 16 21:01:02 PST 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_planarReflectingOperator__
#define __Spheral_planarReflectingOperator__

#include "Geometry/GeomPlane.hh"

namespace Spheral {

template<typename Dimension>
inline
typename Dimension::Tensor
planarReflectingOperator(const typename Dimension::Vector& nhat) {
  return Dimension::Tensor::one() - 2.0*nhat.selfdyad();
}

template<typename Dimension>
inline
typename Dimension::Tensor
planarReflectingOperator(const GeomPlane<Dimension>& plane) {
  return planarReflectingOperator<Dimension>(plane.normal());
}

}

#endif
