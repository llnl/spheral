text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "DEM/SolidBoundary/SphereSolidBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SphereSolidBoundary< Dim< %(ndim)s > >;
}
"""
