text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "DEM/SolidBoundary/InfinitePlaneSolidBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class InfinitePlaneSolidBoundary< Dim< %(ndim)s > >;
}
"""
