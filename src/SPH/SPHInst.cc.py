text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SPH/SPH.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SPH<%(Dim)s>;
}
"""
