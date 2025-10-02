//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "FieldOperations/sampleFieldsMash.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
//============================== sampleFieldsMash() ==============================
template 
FieldList<Dim<1>, Dim<1>::Scalar>sampleFieldsMash<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                         const TableKernel<Dim<1>>& kernel,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& samplePositions,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& sampleWeight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& sampleHfield);

template 
FieldList<Dim<1>, Dim<1>::Vector>sampleFieldsMash<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                         const TableKernel<Dim<1>>& kernel,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& samplePositions,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& sampleWeight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& sampleHfield);

template 
FieldList<Dim<1>, Dim<1>::Tensor>sampleFieldsMash<Dim<1>, Dim<1>::Tensor>(const FieldList<Dim<1>, Dim<1>::Tensor>& fieldList,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                         const TableKernel<Dim<1>>& kernel,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& samplePositions,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& sampleWeight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& sampleHfield);

template 
FieldList<Dim<1>, Dim<1>::SymTensor>sampleFieldsMash<Dim<1>, Dim<1>::SymTensor>(const FieldList<Dim<1>, Dim<1>::SymTensor>& fieldList,
                                            const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                            const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                            const TableKernel<Dim<1>>& kernel,
                                            const FieldList<Dim<1>, Dim<1>::Vector>& samplePositions,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>& sampleWeight,
                                            const FieldList<Dim<1>, Dim<1>::SymTensor>& sampleHfield);
#endif

#if defined(SPHERAL_ENABLE_2D)
//============================== sampleFieldsMash() ==============================
template 
FieldList<Dim<2>, Dim<2>::Scalar>sampleFieldsMash<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                         const TableKernel<Dim<2>>& kernel,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& samplePositions,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& sampleWeight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& sampleHfield);

template 
FieldList<Dim<2>, Dim<2>::Vector>sampleFieldsMash<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                         const TableKernel<Dim<2>>& kernel,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& samplePositions,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& sampleWeight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& sampleHfield);

template 
FieldList<Dim<2>, Dim<2>::Tensor>sampleFieldsMash<Dim<2>, Dim<2>::Tensor>(const FieldList<Dim<2>, Dim<2>::Tensor>& fieldList,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                         const TableKernel<Dim<2>>& kernel,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& samplePositions,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& sampleWeight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& sampleHfield);

template 
FieldList<Dim<2>, Dim<2>::SymTensor>sampleFieldsMash<Dim<2>, Dim<2>::SymTensor>(const FieldList<Dim<2>, Dim<2>::SymTensor>& fieldList,
                                            const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                            const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                            const TableKernel<Dim<2>>& kernel,
                                            const FieldList<Dim<2>, Dim<2>::Vector>& samplePositions,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>& sampleWeight,
                                            const FieldList<Dim<2>, Dim<2>::SymTensor>& sampleHfield);
#endif

#if defined(SPHERAL_ENABLE_3D)
//============================== sampleFieldsMash() ==============================
template 
FieldList<Dim<3>, Dim<3>::Scalar>sampleFieldsMash<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                         const TableKernel<Dim<3>>& kernel,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& samplePositions,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& sampleWeight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& sampleHfield);

template 
FieldList<Dim<3>, Dim<3>::Vector>sampleFieldsMash<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                         const TableKernel<Dim<3>>& kernel,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& samplePositions,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& sampleWeight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& sampleHfield);

template 
FieldList<Dim<3>, Dim<3>::Tensor>sampleFieldsMash<Dim<3>, Dim<3>::Tensor>(const FieldList<Dim<3>, Dim<3>::Tensor>& fieldList,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                         const TableKernel<Dim<3>>& kernel,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& samplePositions,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& sampleWeight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& sampleHfield);

template 
FieldList<Dim<3>, Dim<3>::SymTensor>sampleFieldsMash<Dim<3>, Dim<3>::SymTensor>(const FieldList<Dim<3>, Dim<3>::SymTensor>& fieldList,
                                            const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                            const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                            const TableKernel<Dim<3>>& kernel,
                                            const FieldList<Dim<3>, Dim<3>::Vector>& samplePositions,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>& sampleWeight,
                                            const FieldList<Dim<3>, Dim<3>::SymTensor>& sampleHfield);
#endif

}
