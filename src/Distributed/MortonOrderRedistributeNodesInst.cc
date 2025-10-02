//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Distributed/MortonOrderRedistributeNodes.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class MortonOrderRedistributeNodes<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class MortonOrderRedistributeNodes<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class MortonOrderRedistributeNodes<Dim<3>>;
#endif

}
