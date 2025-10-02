//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Damage/DamageGradientNodeCoupling.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class DamageGradientNodeCoupling<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class DamageGradientNodeCoupling<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class DamageGradientNodeCoupling<Dim<3>>;
#endif

}
