//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/initializeGradients.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template void initializeGradients(const ConnectivityMap<Dim<1>>&, 
                                    const TableKernel<Dim<1>>&, 
                                    const FieldList<Dim<1>, Dim<1>::Vector>&,
                                    const FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                    const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                    const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                    const FieldList<Dim<1>, Dim<1>::Vector>&,
                                          FieldList<Dim<1>, Dim<1>::Tensor>&,
                                          FieldList<Dim<1>, Dim<1>::Vector>&,
                                          FieldList<Dim<1>, Dim<1>::Tensor>&);
#endif

#if defined(SPHERAL_ENABLE_2D)
template void initializeGradients(const ConnectivityMap<Dim<2>>&, 
                                    const TableKernel<Dim<2>>&, 
                                    const FieldList<Dim<2>, Dim<2>::Vector>&,
                                    const FieldList<Dim<2>, Dim<2>::SymTensor>&,
                                    const FieldList<Dim<2>, Dim<2>::Scalar>&,
                                    const FieldList<Dim<2>, Dim<2>::Scalar>&,
                                    const FieldList<Dim<2>, Dim<2>::Vector>&,
                                          FieldList<Dim<2>, Dim<2>::Tensor>&,
                                          FieldList<Dim<2>, Dim<2>::Vector>&,
                                          FieldList<Dim<2>, Dim<2>::Tensor>&);
#endif

#if defined(SPHERAL_ENABLE_3D)
template void initializeGradients(const ConnectivityMap<Dim<3>>&, 
                                    const TableKernel<Dim<3>>&, 
                                    const FieldList<Dim<3>, Dim<3>::Vector>&,
                                    const FieldList<Dim<3>, Dim<3>::SymTensor>&,
                                    const FieldList<Dim<3>, Dim<3>::Scalar>&,
                                    const FieldList<Dim<3>, Dim<3>::Scalar>&,
                                    const FieldList<Dim<3>, Dim<3>::Vector>&,
                                          FieldList<Dim<3>, Dim<3>::Tensor>&,
                                          FieldList<Dim<3>, Dim<3>::Vector>&,
                                          FieldList<Dim<3>, Dim<3>::Tensor>&);
#endif

}
