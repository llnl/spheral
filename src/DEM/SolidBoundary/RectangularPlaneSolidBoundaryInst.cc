//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "DEM/SolidBoundary/RectangularPlaneSolidBoundary.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class RectangularPlaneSolidBoundary<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class RectangularPlaneSolidBoundary<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class RectangularPlaneSolidBoundary<Dim<3>>;
#endif

}
