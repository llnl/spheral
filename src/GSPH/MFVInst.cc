//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/MFV.cc"
#include "GSPH/MFVEvaluateDerivatives.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class MFV<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class MFV<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class MFV<Dim<3>>;
#endif

}
