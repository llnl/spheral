text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ArtificialViscosity/ArtificialViscosityHandle.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ArtificialViscosityHandle<Dim<%(ndim)s>>;
}
"""
