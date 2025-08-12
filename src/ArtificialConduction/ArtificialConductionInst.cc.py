text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ArtificialConduction/ArtificialConduction.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ArtificialConduction< Dim< %(ndim)s > >;
}
"""
