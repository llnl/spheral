//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Damage/JohnsonCookDamage.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class JohnsonCookDamage<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class JohnsonCookDamage<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class JohnsonCookDamage<Dim<3>>;
#endif

}
