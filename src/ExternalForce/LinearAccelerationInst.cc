//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "ExternalForce/LinearAcceleration.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class LinearAcceleration<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class LinearAcceleration<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class LinearAcceleration<Dim<3>>;
#endif

}
