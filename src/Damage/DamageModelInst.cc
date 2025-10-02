//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Damage/DamageModel.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class DamageModel<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class DamageModel<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class DamageModel<Dim<3>>;
#endif

}
