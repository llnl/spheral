
//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "ArtificialViscosity/TensorCRKSPHViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class TensorCRKSPHViscosity<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class TensorCRKSPHViscosity<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class TensorCRKSPHViscosity<Dim<3>>;
#endif

}
