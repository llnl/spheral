//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Boundary/CRKSPHVoidBoundary.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class CRKSPHVoidBoundary<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class CRKSPHVoidBoundary<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class CRKSPHVoidBoundary<Dim<3>>;
#endif

}
