text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SmoothingScale/ASPHClassicSmoothingScale.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ASPHClassicSmoothingScale<Dim<%(ndim)s>>;
}
"""
