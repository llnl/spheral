//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "FieldOperations/sampleMultipleFields2LatticeMash.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
//============================== sampleFieldsMash() ==============================
template 
void
sampleMultipleFields2LatticeMash<Dim<1>>(const FieldListSet<Dim<1>>& fieldListSet,
                                           const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                           const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                           const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                           const FieldList<Dim<1>, int>& mask,
                                           const TableKernel<Dim<1>>& W,
                                           const Dim<1>::Vector& xmin,
                                           const Dim<1>::Vector& xmax,
                                           const vector<int>& nsample,
                                           vector<vector<Dim<1>::Scalar>>& scalarValues,
                                           vector<vector<Dim<1>::Vector>>& vectorValues,
                                           vector<vector<Dim<1>::Tensor>>& tensorValues,
                                           vector<vector<Dim<1>::SymTensor>>& symTensorValues);
#endif

#if defined(SPHERAL_ENABLE_2D)
//============================== sampleFieldsMash() ==============================
template 
void
sampleMultipleFields2LatticeMash<Dim<2>>(const FieldListSet<Dim<2>>& fieldListSet,
                                           const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                           const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                           const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                           const FieldList<Dim<2>, int>& mask,
                                           const TableKernel<Dim<2>>& W,
                                           const Dim<2>::Vector& xmin,
                                           const Dim<2>::Vector& xmax,
                                           const vector<int>& nsample,
                                           vector<vector<Dim<2>::Scalar>>& scalarValues,
                                           vector<vector<Dim<2>::Vector>>& vectorValues,
                                           vector<vector<Dim<2>::Tensor>>& tensorValues,
                                           vector<vector<Dim<2>::SymTensor>>& symTensorValues);
#endif

#if defined(SPHERAL_ENABLE_3D)
//============================== sampleFieldsMash() ==============================
template 
void
sampleMultipleFields2LatticeMash<Dim<3>>(const FieldListSet<Dim<3>>& fieldListSet,
                                           const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                           const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                           const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                           const FieldList<Dim<3>, int>& mask,
                                           const TableKernel<Dim<3>>& W,
                                           const Dim<3>::Vector& xmin,
                                           const Dim<3>::Vector& xmax,
                                           const vector<int>& nsample,
                                           vector<vector<Dim<3>::Scalar>>& scalarValues,
                                           vector<vector<Dim<3>::Vector>>& vectorValues,
                                           vector<vector<Dim<3>::Tensor>>& tensorValues,
                                           vector<vector<Dim<3>::SymTensor>>& symTensorValues);
#endif

}
