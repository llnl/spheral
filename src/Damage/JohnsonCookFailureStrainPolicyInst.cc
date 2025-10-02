//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Damage/JohnsonCookFailureStrainPolicy.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class Spheral::JohnsonCookFailureStrainPolicy<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class Spheral::JohnsonCookFailureStrainPolicy<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class Spheral::JohnsonCookFailureStrainPolicy<Dim<3>>;
#endif

}
