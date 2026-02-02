text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SPH/SolidSPH.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SolidSPH<Dim<%(ndim)s>>;
}
"""
