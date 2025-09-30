
//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "ArtificialViscosity/FiniteVolumeViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class FiniteVolumeViscosity<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class FiniteVolumeViscosity<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class FiniteVolumeViscosity<Dim<3>>;
#endif

}
