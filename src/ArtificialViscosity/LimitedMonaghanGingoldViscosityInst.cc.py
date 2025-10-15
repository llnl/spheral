text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/LimitedMonaghanGingoldViscosity.cc"

namespace Spheral {
  template class LimitedMonaghanGingoldViscosityView< Dim< %(ndim)s > >;
  template class LimitedMonaghanGingoldViscosity< Dim< %(ndim)s > >;
}
"""
