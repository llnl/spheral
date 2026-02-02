text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/Limiters/BarthJespersenLimiter.cc"

namespace Spheral {
  template class BarthJespersenLimiter<Dim< %(ndim)s > >;
}
"""
