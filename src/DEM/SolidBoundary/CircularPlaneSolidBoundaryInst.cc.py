text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "DEM/SolidBoundary/CircularPlaneSolidBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CircularPlaneSolidBoundary< Dim< %(ndim)s > >;
}
"""
