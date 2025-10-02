//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/RiemannSolvers/SecondOrderArtificialViscosity.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class SecondOrderArtificialViscosity<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class SecondOrderArtificialViscosity<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class SecondOrderArtificialViscosity<Dim<3>>;
#endif

}
