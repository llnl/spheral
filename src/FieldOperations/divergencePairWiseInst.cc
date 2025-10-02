//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "FieldOperations/divergencePairWise.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
//============================== gradient() ==============================
template 
FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Vector>::DivergenceType>divergencePairWise<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                           const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                           const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                           const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                           const FieldList<Dim<1>, Dim<1>::Scalar>& density,
                                           const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                           const TableKernel<Dim<1>>& kernel);
#endif

#if defined(SPHERAL_ENABLE_2D)
//============================== gradient() ==============================
template 
FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Vector>::DivergenceType>divergencePairWise<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                           const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                           const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                           const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                           const FieldList<Dim<2>, Dim<2>::Scalar>& density,
                                           const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                           const TableKernel<Dim<2>>& kernel);
#endif

#if defined(SPHERAL_ENABLE_3D)
//============================== gradient() ==============================
template 
FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Vector>::DivergenceType>divergencePairWise<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                           const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                           const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                           const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                           const FieldList<Dim<3>, Dim<3>::Scalar>& density,
                                           const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                           const TableKernel<Dim<3>>& kernel);
#endif

}
