//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Distributed/NestedGridDistributedBoundary.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class NestedGridDistributedBoundary<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class NestedGridDistributedBoundary<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class NestedGridDistributedBoundary<Dim<3>>;
#endif

}
