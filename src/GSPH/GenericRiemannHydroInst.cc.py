text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GSPH/GenericRiemannHydro.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class GenericRiemannHydro< Dim< %(ndim)s > >;
}
"""
