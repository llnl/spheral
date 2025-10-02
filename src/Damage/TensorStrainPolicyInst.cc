//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Damage/TensorStrainPolicy.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class Spheral::TensorStrainPolicy<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class Spheral::TensorStrainPolicy<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class Spheral::TensorStrainPolicy<Dim<3>>;
#endif

}
