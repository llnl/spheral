//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "CRKSPH/CRKSPHBase.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class CRKSPHBase<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class CRKSPHBase<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class CRKSPHBase<Dim<3>>;
#endif

}
