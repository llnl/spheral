//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Distributed/SpaceFillingCurveRedistributeNodes.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class SpaceFillingCurveRedistributeNodes<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class SpaceFillingCurveRedistributeNodes<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class SpaceFillingCurveRedistributeNodes<Dim<3>>;
#endif

}
