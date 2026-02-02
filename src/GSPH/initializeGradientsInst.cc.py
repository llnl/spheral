text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GSPH/initializeGradients.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {


  template void initializeGradients(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                    const TableKernel<Dim< %(ndim)s > >&, 
                                    const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                    const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                    const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                    const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                    const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                          FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                          FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                          FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&);
}
"""
