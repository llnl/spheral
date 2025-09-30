
//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "ArtificialViscosity/MonaghanGingoldViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class MonaghanGingoldViscosity<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class MonaghanGingoldViscosity<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class MonaghanGingoldViscosity<Dim<3>>;
#endif

}
