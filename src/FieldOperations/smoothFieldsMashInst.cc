//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "FieldOperations/smoothFieldsMash.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template 
FieldList<Dim<1>, Dim<1>::Scalar>smoothFieldsMash<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                         const TableKernel<Dim<1>>& kernel);
template 
FieldList<Dim<1>, Dim<1>::Vector>smoothFieldsMash<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                         const TableKernel<Dim<1>>& kernel);
template 
FieldList<Dim<1>, Dim<1>::Tensor>smoothFieldsMash<Dim<1>, Dim<1>::Tensor>(const FieldList<Dim<1>, Dim<1>::Tensor>& fieldList,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                         const TableKernel<Dim<1>>& kernel);
template 
FieldList<Dim<1>, Dim<1>::SymTensor>smoothFieldsMash<Dim<1>, Dim<1>::SymTensor>(const FieldList<Dim<1>, Dim<1>::SymTensor>& fieldList,
                                            const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                            const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                            const TableKernel<Dim<1>>& kernel);
#endif

#if defined(SPHERAL_ENABLE_2D)
template 
FieldList<Dim<2>, Dim<2>::Scalar>smoothFieldsMash<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                         const TableKernel<Dim<2>>& kernel);
template 
FieldList<Dim<2>, Dim<2>::Vector>smoothFieldsMash<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                         const TableKernel<Dim<2>>& kernel);
template 
FieldList<Dim<2>, Dim<2>::Tensor>smoothFieldsMash<Dim<2>, Dim<2>::Tensor>(const FieldList<Dim<2>, Dim<2>::Tensor>& fieldList,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                         const TableKernel<Dim<2>>& kernel);
template 
FieldList<Dim<2>, Dim<2>::SymTensor>smoothFieldsMash<Dim<2>, Dim<2>::SymTensor>(const FieldList<Dim<2>, Dim<2>::SymTensor>& fieldList,
                                            const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                            const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                            const TableKernel<Dim<2>>& kernel);
#endif

#if defined(SPHERAL_ENABLE_3D)
template 
FieldList<Dim<3>, Dim<3>::Scalar>smoothFieldsMash<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                         const TableKernel<Dim<3>>& kernel);
template 
FieldList<Dim<3>, Dim<3>::Vector>smoothFieldsMash<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                         const TableKernel<Dim<3>>& kernel);
template 
FieldList<Dim<3>, Dim<3>::Tensor>smoothFieldsMash<Dim<3>, Dim<3>::Tensor>(const FieldList<Dim<3>, Dim<3>::Tensor>& fieldList,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                         const TableKernel<Dim<3>>& kernel);
template 
FieldList<Dim<3>, Dim<3>::SymTensor>smoothFieldsMash<Dim<3>, Dim<3>::SymTensor>(const FieldList<Dim<3>, Dim<3>::SymTensor>& fieldList,
                                            const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                            const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                            const TableKernel<Dim<3>>& kernel);
#endif

}
