//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/TensorMonaghanGingoldViscosity.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class TensorMonaghanGingoldViscosity<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class TensorMonaghanGingoldViscosity<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class TensorMonaghanGingoldViscosity<Dim<3>>;
#endif

}
