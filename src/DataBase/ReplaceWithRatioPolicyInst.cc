//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "DataBase/ReplaceWithRatioPolicy.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class ReplaceWithRatioPolicy<Dim<1>, Dim<1>::Scalar>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class ReplaceWithRatioPolicy<Dim<2>, Dim<2>::Scalar>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class ReplaceWithRatioPolicy<Dim<3>, Dim<3>::Scalar>;
#endif

}
