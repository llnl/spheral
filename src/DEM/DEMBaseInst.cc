//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "DEM/DEMBase.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class DEMBase<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class DEMBase<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class DEMBase<Dim<3>>;
#endif

}
