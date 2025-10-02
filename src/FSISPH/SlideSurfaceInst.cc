//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "FSISPH/SlideSurface.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class SlideSurface<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class SlideSurface<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class SlideSurface<Dim<3>>;
#endif

}
