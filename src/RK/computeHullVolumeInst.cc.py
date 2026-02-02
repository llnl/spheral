text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "RK/computeHullVolume.cc"

namespace Spheral {
  template void computeHullVolume(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                  const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                  const ConnectivityMap<Dim<%(ndim)s>>&, 
                                  const bool,
                                  FieldList<Dim< %(ndim)s >, int>&,
                                  FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                  FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::FacetedVolume>&);
}

"""
