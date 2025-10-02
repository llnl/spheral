//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "DEM/LinearSpringDEM.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class LinearSpringDEM<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class LinearSpringDEM<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class LinearSpringDEM<Dim<3>>;
#endif

}
