text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SmoothingScale/SmoothingScaleBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SmoothingScaleBase<Dim<%(ndim)s>>;
}
"""
