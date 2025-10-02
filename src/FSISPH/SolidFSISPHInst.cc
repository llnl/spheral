//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "FSISPH/SolidFSISPH.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class SolidFSISPH<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class SolidFSISPH<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class SolidFSISPH<Dim<3>>;
#endif

}
