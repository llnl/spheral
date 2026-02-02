text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Integrator/ForwardEuler.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ForwardEuler< Dim< %(ndim)s > >;
}
"""
