//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Damage/StrainPolicy.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class Spheral::StrainPolicy<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class Spheral::StrainPolicy<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class Spheral::StrainPolicy<Dim<3>>;
#endif

}
