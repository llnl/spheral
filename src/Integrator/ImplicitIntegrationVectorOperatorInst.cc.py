text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Integrator/ImplicitIntegrationVectorOperator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ImplicitIntegrationVectorOperator<Dim<%(ndim)s>>;
}

"""
