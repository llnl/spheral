text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/SingularityEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SingularityEquationOfState< Dim< %(ndim)s > >;
}
"""
