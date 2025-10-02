//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Distributed/DistributeByXPosition.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class DistributeByXPosition<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class DistributeByXPosition<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class DistributeByXPosition<Dim<3>>;
#endif

}
