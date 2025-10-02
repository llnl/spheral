//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "FieldOperations/gradient.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
//============================== gradient() ==============================
template 
FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Scalar>::GradientType>gradient<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                 const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                 const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                 const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                 const FieldList<Dim<1>, Dim<1>::Scalar>& density,
                                 const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                 const TableKernel<Dim<1>>& kernel);
template 
FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Vector>::GradientType>gradient<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                 const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                 const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                 const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                 const FieldList<Dim<1>, Dim<1>::Scalar>& density,
                                 const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                 const TableKernel<Dim<1>>& kernel);

template 
FieldList<Dim<1>, std::vector<MathTraits<Dim<1>, Dim<1>::Scalar>::GradientType>>gradient<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, std::vector<Dim<1>::Scalar>>& fieldList,
                                 const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                 const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                 const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                 const FieldList<Dim<1>, Dim<1>::Scalar>& density,
                                 const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                 const TableKernel<Dim<1>>& kernel);
template 
FieldList<Dim<1>, std::vector<MathTraits<Dim<1>, Dim<1>::Vector>::GradientType>>gradient<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, std::vector<Dim<1>::Vector>>& fieldList,
                                 const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                 const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                 const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                 const FieldList<Dim<1>, Dim<1>::Scalar>& density,
                                 const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                 const TableKernel<Dim<1>>& kernel);


//============================== limiter() ==============================
template 
FieldList<Dim<1>, Dim<1>::SymTensor>limiter<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                const FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Scalar>::GradientType>& gradient,
                                const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                const TableKernel<Dim<1>>& kernel);
template 
FieldList<Dim<1>, Dim<1>::SymTensor>limiter<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                const FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Vector>::GradientType>& gradient,
                                const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                const TableKernel<Dim<1>>& kernel);
#endif

#if defined(SPHERAL_ENABLE_2D)
//============================== gradient() ==============================
template 
FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Scalar>::GradientType>gradient<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                                 const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                 const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                 const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                 const FieldList<Dim<2>, Dim<2>::Scalar>& density,
                                 const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                 const TableKernel<Dim<2>>& kernel);
template 
FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Vector>::GradientType>gradient<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                 const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                 const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                 const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                 const FieldList<Dim<2>, Dim<2>::Scalar>& density,
                                 const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                 const TableKernel<Dim<2>>& kernel);

template 
FieldList<Dim<2>, std::vector<MathTraits<Dim<2>, Dim<2>::Scalar>::GradientType>>gradient<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, std::vector<Dim<2>::Scalar>>& fieldList,
                                 const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                 const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                 const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                 const FieldList<Dim<2>, Dim<2>::Scalar>& density,
                                 const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                 const TableKernel<Dim<2>>& kernel);
template 
FieldList<Dim<2>, std::vector<MathTraits<Dim<2>, Dim<2>::Vector>::GradientType>>gradient<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, std::vector<Dim<2>::Vector>>& fieldList,
                                 const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                 const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                 const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                 const FieldList<Dim<2>, Dim<2>::Scalar>& density,
                                 const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                 const TableKernel<Dim<2>>& kernel);


//============================== limiter() ==============================
template 
FieldList<Dim<2>, Dim<2>::SymTensor>limiter<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                                const FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Scalar>::GradientType>& gradient,
                                const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                const TableKernel<Dim<2>>& kernel);
template 
FieldList<Dim<2>, Dim<2>::SymTensor>limiter<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                const FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Vector>::GradientType>& gradient,
                                const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                const TableKernel<Dim<2>>& kernel);
#endif

#if defined(SPHERAL_ENABLE_3D)
//============================== gradient() ==============================
template 
FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Scalar>::GradientType>gradient<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                                 const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                 const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                 const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                 const FieldList<Dim<3>, Dim<3>::Scalar>& density,
                                 const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                 const TableKernel<Dim<3>>& kernel);
template 
FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Vector>::GradientType>gradient<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                 const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                 const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                 const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                 const FieldList<Dim<3>, Dim<3>::Scalar>& density,
                                 const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                 const TableKernel<Dim<3>>& kernel);

template 
FieldList<Dim<3>, std::vector<MathTraits<Dim<3>, Dim<3>::Scalar>::GradientType>>gradient<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, std::vector<Dim<3>::Scalar>>& fieldList,
                                 const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                 const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                 const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                 const FieldList<Dim<3>, Dim<3>::Scalar>& density,
                                 const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                 const TableKernel<Dim<3>>& kernel);
template 
FieldList<Dim<3>, std::vector<MathTraits<Dim<3>, Dim<3>::Vector>::GradientType>>gradient<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, std::vector<Dim<3>::Vector>>& fieldList,
                                 const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                 const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                 const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                 const FieldList<Dim<3>, Dim<3>::Scalar>& density,
                                 const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                 const TableKernel<Dim<3>>& kernel);


//============================== limiter() ==============================
template 
FieldList<Dim<3>, Dim<3>::SymTensor>limiter<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                                const FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Scalar>::GradientType>& gradient,
                                const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                const TableKernel<Dim<3>>& kernel);
template 
FieldList<Dim<3>, Dim<3>::SymTensor>limiter<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                const FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Vector>::GradientType>& gradient,
                                const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                const TableKernel<Dim<3>>& kernel);
#endif

}
