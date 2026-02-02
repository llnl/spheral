text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Integrator/CrankNicolson.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CrankNicolson<Dim<%(ndim)s>>;
}
"""
