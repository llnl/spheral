text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/MonaghanGingoldViscosity.cc"

namespace Spheral {
  template class MonaghanGingoldViscosityView< Dim< %(ndim)s > >;
  template class MonaghanGingoldViscosity< Dim< %(ndim)s > >;
}
"""
