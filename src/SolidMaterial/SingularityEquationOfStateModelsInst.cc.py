text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/SingularityEquationOfStateModels.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SingularityGammaLawGas< Dim< %(ndim)s > >;
  template class SingularityExtendedVinet< Dim< %(ndim)s > >;
  template class SingularitySpiner< Dim< %(ndim)s > >;
  template class SingularitySpinerT< Dim< %(ndim)s > >;
  template class SingularityEOSPAC< Dim< %(ndim)s > >;
}
"""
