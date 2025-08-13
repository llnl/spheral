text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SPH/PSPH.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PSPH<Dim<%(ndim)s>>;
}
"""
