//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Damage/TensorDamagePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class TensorDamagePolicy<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class TensorDamagePolicy<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class TensorDamagePolicy<Dim<3>>;
#endif

}
