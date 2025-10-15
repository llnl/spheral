text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/MonaghanGingoldViscosity.cc"

namespace Spheral {
  template class MonaghanGingoldViscosity< Dim< %(ndim)s > >;
}
"""
