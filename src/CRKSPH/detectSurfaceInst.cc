//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "CRKSPH/detectSurface.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template void detectSurface(const ConnectivityMap<Dim<1>>& connectivityMap,
                            const FieldList<Dim<1>, Dim<1>::Scalar>& m0,
                            const FieldList<Dim<1>, Dim<1>::Vector>& m1,
                            const FieldList<Dim<1>, Dim<1>::Vector>& position,
                            const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                            const double detectThreshold,
                            const double detectRange,
                            const double sweepAngle,
                            FieldList<Dim<1>, int>& surfacePoint);
#endif

#if defined(SPHERAL_ENABLE_2D)
template void detectSurface(const ConnectivityMap<Dim<2>>& connectivityMap,
                            const FieldList<Dim<2>, Dim<2>::Scalar>& m0,
                            const FieldList<Dim<2>, Dim<2>::Vector>& m1,
                            const FieldList<Dim<2>, Dim<2>::Vector>& position,
                            const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                            const double detectThreshold,
                            const double detectRange,
                            const double sweepAngle,
                            FieldList<Dim<2>, int>& surfacePoint);
#endif

#if defined(SPHERAL_ENABLE_3D)
template void detectSurface(const ConnectivityMap<Dim<3>>& connectivityMap,
                            const FieldList<Dim<3>, Dim<3>::Scalar>& m0,
                            const FieldList<Dim<3>, Dim<3>::Vector>& m1,
                            const FieldList<Dim<3>, Dim<3>::Vector>& position,
                            const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                            const double detectThreshold,
                            const double detectRange,
                            const double sweepAngle,
                            FieldList<Dim<3>, int>& surfacePoint);
#endif

}
