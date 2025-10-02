//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "ExternalForce/ConstantAcceleration.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class ConstantAcceleration<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class ConstantAcceleration<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class ConstantAcceleration<Dim<3>>;
#endif

}
