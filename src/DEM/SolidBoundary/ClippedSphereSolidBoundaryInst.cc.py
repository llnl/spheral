text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "DEM/SolidBoundary/ClippedSphereSolidBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ClippedSphereSolidBoundary< Dim< %(ndim)s > >;
}
"""
