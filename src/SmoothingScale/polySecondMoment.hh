//------------------------------------------------------------------------------
// Compute the second moment about the give position for a polytope
//
// Note, these methods currently assume the polytopes are convex.
//------------------------------------------------------------------------------
#ifndef __Spheral_polySecondMoment__
#define __Spheral_polySecondMoment__

#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// 1D
//------------------------------------------------------------------------------
Dim<1>::SymTensor
polySecondMoment(const Dim<1>::FacetedVolume& poly,
                 const Dim<1>::Vector& center);

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
Dim<2>::SymTensor
polySecondMoment(const Dim<2>::FacetedVolume& poly,
                 const Dim<2>::Vector& center);

//------------------------------------------------------------------------------
// 3D
//------------------------------------------------------------------------------
Dim<3>::SymTensor
polySecondMoment(const Dim<3>::FacetedVolume& poly,
                 const Dim<3>::Vector& center);

}

#endif
