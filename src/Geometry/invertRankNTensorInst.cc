//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/invertRankNTensor.cc"
#include "Geometry/GeomFourthRankTensor.hh"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template GeomFourthRankTensor<1> invertRankNTensor(const GeomFourthRankTensor<1>&);
#endif

#if defined(SPHERAL_ENABLE_2D)
template GeomFourthRankTensor<2> invertRankNTensor(const GeomFourthRankTensor<2>&);
#endif

#if defined(SPHERAL_ENABLE_3D)
template GeomFourthRankTensor<3> invertRankNTensor(const GeomFourthRankTensor<3>&);
#endif

}
