//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/Limiters/OspreLimiter.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class OspreLimiter<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class OspreLimiter<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class OspreLimiter<Dim<3>>;
#endif

}
