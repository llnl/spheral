text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SmoothingScale/SPHSmoothingScale.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SPHSmoothingScale<Dim<%(ndim)s>>;
}
"""
