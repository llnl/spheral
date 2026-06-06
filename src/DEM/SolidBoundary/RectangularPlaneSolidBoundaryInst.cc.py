text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "DEM/SolidBoundary/RectangularPlaneSolidBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class RectangularPlaneSolidBoundary< Dim< %(ndim)s > >;
}
"""
