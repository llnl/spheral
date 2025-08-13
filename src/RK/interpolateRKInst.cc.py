text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "RK/interpolateRK.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

template 
std::vector<std::variant<FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>,
                         FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>,
                         FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>,
                         FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::SymTensor>,
                         FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::ThirdRankTensor>>>
interpolateRK<Dim<%(ndim)s>>(const std::vector<std::variant<FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>,
                                                            FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>,
                                                            FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>,
                                                            FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::SymTensor>,
                                                            FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::ThirdRankTensor>>>& fieldLists,
                                          const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>& position,
                                          const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>& weight,
                                          const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::SymTensor>& H,
                                          const ConnectivityMap<Dim<%(ndim)s> >& connectivityMap,
                                          const ReproducingKernel< Dim<%(ndim)s> >& kernel,
                                          const FieldList<Dim<%(ndim)s>, RKCoefficients<Dim<%(ndim)s>>>& corrections,
                                          const NodeCoupling& nodeCoupling);

}

"""
