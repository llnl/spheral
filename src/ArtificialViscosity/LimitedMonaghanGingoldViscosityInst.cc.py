text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/LimitedMonaghanGingoldViscosityView.cc"

namespace Spheral {
  template class LimitedMonaghanGingoldViscosityView< Dim< %(ndim)s > >;
}
"""
