text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SmoothingScale/IncrementASPHHtensor.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class IncrementASPHHtensor<Dim<%(ndim)s>>;
}
"""
