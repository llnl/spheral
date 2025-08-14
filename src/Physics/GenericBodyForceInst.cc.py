text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Physics/GenericBodyForce.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class GenericBodyForce< Dim< %(ndim)s > >;
}
"""
