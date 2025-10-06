text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FSISPH/computeFSISPHSumMassDensity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {


  template void computeFSISPHSumMassDensity(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                            const TableKernel<Dim< %(ndim)s > >&, 
                                            const std::vector<int>&,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                            const bool consistentSum,
                                                  FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
}
"""
