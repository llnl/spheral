text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GSPH/Policies/MassFluxPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MassFluxPolicy<Dim< %(ndim)s > >;
}
"""
