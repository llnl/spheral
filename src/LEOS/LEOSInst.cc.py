text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "LEOS.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class LEOS<Dim< %(ndim)s > >;
}
"""
