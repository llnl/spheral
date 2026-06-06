//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------
#include "PolyGravity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

#ifdef SPHERAL2D
  template class PolyGravity<Dim<2> >;
#endif
#ifdef SPHERAL3D
  template class PolyGravity<Dim<3> >;
#endif

}
