text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/LimitedMonaghanGingoldViscosity.hh"

namespace Spheral {
  template class LimitedMonaghanGingoldViscosity< Dim< %(ndim)s > >;
  template class LimitedMonaghanGingoldViscosityView< Dim< %(ndim)s > >;
}
"""
