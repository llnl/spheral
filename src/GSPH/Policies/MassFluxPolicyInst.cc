//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/Policies/MassFluxPolicy.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class MassFluxPolicy<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class MassFluxPolicy<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class MassFluxPolicy<Dim<3>>;
#endif

}
