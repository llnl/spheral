text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/MonaghanGingoldViscosity.hh"

namespace Spheral {
  template class MonaghanGingoldViscosity< Dim< %(ndim)s > >;
  template class MonaghanGingoldViscosityView< Dim< %(ndim)s > >;
}
"""
