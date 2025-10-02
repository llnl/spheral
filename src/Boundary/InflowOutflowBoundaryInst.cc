//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Boundary/InflowOutflowBoundary.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class InflowOutflowBoundary<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class InflowOutflowBoundary<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class InflowOutflowBoundary<Dim<3>>;
#endif

}
