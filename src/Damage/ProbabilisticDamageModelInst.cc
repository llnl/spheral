//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Damage/ProbabilisticDamageModel.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class ProbabilisticDamageModel<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class ProbabilisticDamageModel<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class ProbabilisticDamageModel<Dim<3>>;
#endif

}
