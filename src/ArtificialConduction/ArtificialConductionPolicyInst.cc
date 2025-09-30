
//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "ArtificialConduction/ArtificialConductionPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class ArtificialConductionPolicy<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class ArtificialConductionPolicy<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class ArtificialConductionPolicy<Dim<3>>;
#endif

}
