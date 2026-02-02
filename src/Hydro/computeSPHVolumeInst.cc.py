text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Hydro/computeSPHVolume.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {


  template void computeSPHVolume(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                 const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                       FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
}
"""
