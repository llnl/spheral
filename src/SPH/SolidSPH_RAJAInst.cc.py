text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SPH/SolidSPH_RAJA.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SolidSPH_RAJA<Dim<%(ndim)s>>;
}
"""
