//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Boundary/ReflectingBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class ReflectingBoundary<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class ReflectingBoundary<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class ReflectingBoundary<Dim<3>>;
#endif

}
