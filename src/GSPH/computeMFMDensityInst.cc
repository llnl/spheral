//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/computeMFMDensity.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template void computeMFMDensity(const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                  const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                        FieldList<Dim<1>, Dim<1>::Scalar>&);
#endif

#if defined(SPHERAL_ENABLE_2D)
template void computeMFMDensity(const FieldList<Dim<2>, Dim<2>::Scalar>&,
                                  const FieldList<Dim<2>, Dim<2>::Scalar>&,
                                        FieldList<Dim<2>, Dim<2>::Scalar>&);
#endif

#if defined(SPHERAL_ENABLE_3D)
template void computeMFMDensity(const FieldList<Dim<3>, Dim<3>::Scalar>&,
                                  const FieldList<Dim<3>, Dim<3>::Scalar>&,
                                        FieldList<Dim<3>, Dim<3>::Scalar>&);
#endif

}
