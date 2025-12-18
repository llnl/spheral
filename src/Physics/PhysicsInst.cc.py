text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Physics/Physics.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class Physics<Dim< %(ndim)s > >;
}
"""
