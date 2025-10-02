//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Boundary/RigidBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class RigidBoundary<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class RigidBoundary<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class RigidBoundary<Dim<3>>;
#endif

}
