text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SPH/SPHBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SPHBase<%(Dim)s>;
}
"""
