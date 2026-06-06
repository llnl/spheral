text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Integrator/ImplicitIntegrator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ImplicitIntegrator<Dim<%(ndim)s>>;
}
"""
