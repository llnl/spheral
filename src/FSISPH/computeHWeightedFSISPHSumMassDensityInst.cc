//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "FSISPH/computeHWeightedFSISPHSumMassDensity.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template void computeHWeightedFSISPHSumMassDensity(const ConnectivityMap<Dim<1>>&, 
                                                     const TableKernel<Dim<1>>&, 
                                                     const std::vector<int>&,
                                                     const FieldList<Dim<1>, Dim<1>::Vector>&,
                                                     const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                                     const FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                                           FieldList<Dim<1>, Dim<1>::Scalar>&);
#endif

#if defined(SPHERAL_ENABLE_2D)
template void computeHWeightedFSISPHSumMassDensity(const ConnectivityMap<Dim<2>>&, 
                                                     const TableKernel<Dim<2>>&, 
                                                     const std::vector<int>&,
                                                     const FieldList<Dim<2>, Dim<2>::Vector>&,
                                                     const FieldList<Dim<2>, Dim<2>::Scalar>&,
                                                     const FieldList<Dim<2>, Dim<2>::SymTensor>&,
                                                           FieldList<Dim<2>, Dim<2>::Scalar>&);
#endif

#if defined(SPHERAL_ENABLE_3D)
template void computeHWeightedFSISPHSumMassDensity(const ConnectivityMap<Dim<3>>&, 
                                                     const TableKernel<Dim<3>>&, 
                                                     const std::vector<int>&,
                                                     const FieldList<Dim<3>, Dim<3>::Vector>&,
                                                     const FieldList<Dim<3>, Dim<3>::Scalar>&,
                                                     const FieldList<Dim<3>, Dim<3>::SymTensor>&,
                                                           FieldList<Dim<3>, Dim<3>::Scalar>&);
#endif

}
