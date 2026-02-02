text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/RiemannSolvers/SecondOrderArtificialViscosity.cc"

namespace Spheral {
  template class SecondOrderArtificialViscosity<Dim< %(ndim)s > >;
}
"""
