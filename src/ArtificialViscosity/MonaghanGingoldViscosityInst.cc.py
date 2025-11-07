text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/MonaghanGingoldViscosityView.cc"

namespace Spheral {
  template class MonaghanGingoldViscosityView< Dim< %(ndim)s > >;
}
"""
