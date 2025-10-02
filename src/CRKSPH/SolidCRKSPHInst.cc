//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "CRKSPH/SolidCRKSPH.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class SolidCRKSPH<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class SolidCRKSPH<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class SolidCRKSPH<Dim<3>>;
#endif

}
