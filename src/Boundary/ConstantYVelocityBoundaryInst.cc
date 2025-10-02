//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Boundary/ConstantYVelocityBoundary.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class ConstantYVelocityBoundary<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class ConstantYVelocityBoundary<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class ConstantYVelocityBoundary<Dim<3>>;
#endif

}
