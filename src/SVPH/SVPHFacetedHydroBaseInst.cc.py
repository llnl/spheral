text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SVPH/SVPHFacetedHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SVPHFacetedHydroBase< Dim< %(ndim)s > >;
}
"""
