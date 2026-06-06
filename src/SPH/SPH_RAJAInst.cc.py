text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SPH/SPH_RAJA.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SPH_RAJA<%(Dim)s>;
}
"""
