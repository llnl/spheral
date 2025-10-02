//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Gravity/NBodyGravity.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class NBodyGravity<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class NBodyGravity<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class NBodyGravity<Dim<3>>;
#endif

}
