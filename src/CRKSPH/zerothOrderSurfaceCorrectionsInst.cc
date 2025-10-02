//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "CRKSPH/zerothOrderSurfaceCorrections.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template void
zerothOrderSurfaceCorrections(FieldList<Dim<1>, Dim<1>::Scalar>& A,
                              FieldList<Dim<1>, Dim<1>::Vector>& B,
                              FieldList<Dim<1>, Dim<1>::Tensor>& C,
                              FieldList<Dim<1>, Dim<1>::Vector>& gradA,
                              FieldList<Dim<1>, Dim<1>::Tensor>& gradB,
                              FieldList<Dim<1>, Dim<1>::ThirdRankTensor>& gradC,
                              const FieldList<Dim<1>, Dim<1>::Scalar>& m0,
                              const FieldList<Dim<1>, Dim<1>::Vector>& gradm0,
                              const FieldList<Dim<1>, int>& surfacePoint);
#endif

#if defined(SPHERAL_ENABLE_2D)
template void
zerothOrderSurfaceCorrections(FieldList<Dim<2>, Dim<2>::Scalar>& A,
                              FieldList<Dim<2>, Dim<2>::Vector>& B,
                              FieldList<Dim<2>, Dim<2>::Tensor>& C,
                              FieldList<Dim<2>, Dim<2>::Vector>& gradA,
                              FieldList<Dim<2>, Dim<2>::Tensor>& gradB,
                              FieldList<Dim<2>, Dim<2>::ThirdRankTensor>& gradC,
                              const FieldList<Dim<2>, Dim<2>::Scalar>& m0,
                              const FieldList<Dim<2>, Dim<2>::Vector>& gradm0,
                              const FieldList<Dim<2>, int>& surfacePoint);
#endif

#if defined(SPHERAL_ENABLE_3D)
template void
zerothOrderSurfaceCorrections(FieldList<Dim<3>, Dim<3>::Scalar>& A,
                              FieldList<Dim<3>, Dim<3>::Vector>& B,
                              FieldList<Dim<3>, Dim<3>::Tensor>& C,
                              FieldList<Dim<3>, Dim<3>::Vector>& gradA,
                              FieldList<Dim<3>, Dim<3>::Tensor>& gradB,
                              FieldList<Dim<3>, Dim<3>::ThirdRankTensor>& gradC,
                              const FieldList<Dim<3>, Dim<3>::Scalar>& m0,
                              const FieldList<Dim<3>, Dim<3>::Vector>& gradm0,
                              const FieldList<Dim<3>, int>& surfacePoint);
#endif

}
