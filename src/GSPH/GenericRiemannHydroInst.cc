//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/GenericRiemannHydro.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class GenericRiemannHydro<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class GenericRiemannHydro<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class GenericRiemannHydro<Dim<3>>;
#endif

}
