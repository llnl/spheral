text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FSISPH/SlideSurface.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SlideSurface< Dim< %(ndim)s > >;
}
"""
