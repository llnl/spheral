//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/Limiters/VanAlbaLimiter.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class VanAlbaLimiter<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class VanAlbaLimiter<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class VanAlbaLimiter<Dim<3>>;
#endif

}
