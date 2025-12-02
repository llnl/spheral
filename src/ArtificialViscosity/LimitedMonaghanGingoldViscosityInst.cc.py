text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/LimitedMonaghanGingoldViscosityView.hh"

namespace Spheral {
  template class LimitedMonaghanGingoldViscosityView< Dim< %(ndim)s > >;
}
"""
