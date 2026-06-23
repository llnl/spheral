text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ArtificialViscosity/ArtificialViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ArtificialViscosity<Dim<%(ndim)s>>;
}
"""
