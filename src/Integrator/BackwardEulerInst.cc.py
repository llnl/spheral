text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Integrator/BackwardEuler.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class BackwardEuler<Dim<%(ndim)s>>;
}
"""
