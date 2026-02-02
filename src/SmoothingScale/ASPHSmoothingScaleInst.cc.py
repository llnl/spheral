text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SmoothingScale/ASPHSmoothingScale.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ASPHSmoothingScale<Dim<%(ndim)s>>;
}
"""
