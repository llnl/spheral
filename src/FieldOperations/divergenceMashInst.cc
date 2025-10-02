//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "FieldOperations/divergenceMash.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
//============================== divergenceMash() ==============================
template 
FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Vector>::DivergenceType>divergenceMash<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                       const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                       const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                       const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                       const TableKernel<Dim<1>>& kernel);

template 
FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Tensor>::DivergenceType>divergenceMash<Dim<1>, Dim<1>::Tensor>(const FieldList<Dim<1>, Dim<1>::Tensor>& fieldList,
                                       const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                       const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                       const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                       const TableKernel<Dim<1>>& kernel);
#endif

#if defined(SPHERAL_ENABLE_2D)
//============================== divergenceMash() ==============================
template 
FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Vector>::DivergenceType>divergenceMash<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                       const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                       const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                       const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                       const TableKernel<Dim<2>>& kernel);

template 
FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Tensor>::DivergenceType>divergenceMash<Dim<2>, Dim<2>::Tensor>(const FieldList<Dim<2>, Dim<2>::Tensor>& fieldList,
                                       const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                       const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                       const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                       const TableKernel<Dim<2>>& kernel);
#endif

#if defined(SPHERAL_ENABLE_3D)
//============================== divergenceMash() ==============================
template 
FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Vector>::DivergenceType>divergenceMash<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                       const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                       const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                       const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                       const TableKernel<Dim<3>>& kernel);

template 
FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Tensor>::DivergenceType>divergenceMash<Dim<3>, Dim<3>::Tensor>(const FieldList<Dim<3>, Dim<3>::Tensor>& fieldList,
                                       const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                       const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                       const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                       const TableKernel<Dim<3>>& kernel);
#endif

}
