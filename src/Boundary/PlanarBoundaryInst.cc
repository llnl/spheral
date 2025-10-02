//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Boundary/PlanarBoundary.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class PlanarBoundary<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class PlanarBoundary<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class PlanarBoundary<Dim<3>>;
#endif

}
