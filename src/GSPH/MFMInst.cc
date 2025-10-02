//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/MFM.cc"
#include "GSPH/MFMEvaluateDerivatives.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class MFM<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class MFM<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class MFM<Dim<3>>;
#endif

}
