//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/TensorSVPHViscosity.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class TensorSVPHViscosity<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class TensorSVPHViscosity<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class TensorSVPHViscosity<Dim<3>>;
#endif

}
