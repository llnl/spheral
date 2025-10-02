//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/GSPH.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class GSPH<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class GSPH<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class GSPH<Dim<3>>;
#endif

}
