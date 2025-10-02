//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Distributed/VoronoiRedistributeNodes.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class VoronoiRedistributeNodes<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class VoronoiRedistributeNodes<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class VoronoiRedistributeNodes<Dim<3>>;
#endif

}
