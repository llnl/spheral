text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/MonaghanGingoldViscosityView.hh"

namespace Spheral {
  template class MonaghanGingoldViscosityView< Dim< %(ndim)s > >;
}
"""
