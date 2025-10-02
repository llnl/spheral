//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "DEM/SolidBoundary/CylinderSolidBoundary.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class CylinderSolidBoundary<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class CylinderSolidBoundary<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class CylinderSolidBoundary<Dim<3>>;
#endif

}
