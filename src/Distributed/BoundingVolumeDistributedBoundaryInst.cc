//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Distributed/BoundingVolumeDistributedBoundary.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class BoundingVolumeDistributedBoundary<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class BoundingVolumeDistributedBoundary<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class BoundingVolumeDistributedBoundary<Dim<3>>;
#endif

}
