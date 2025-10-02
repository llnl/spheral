//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "DEM/SolidBoundary/SolidBoundaryBase.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class SolidBoundaryBase<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class SolidBoundaryBase<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class SolidBoundaryBase<Dim<3>>;
#endif

}
