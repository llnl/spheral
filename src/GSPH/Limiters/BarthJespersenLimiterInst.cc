//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/Limiters/BarthJespersenLimiter.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class BarthJespersenLimiter<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class BarthJespersenLimiter<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class BarthJespersenLimiter<Dim<3>>;
#endif

}
