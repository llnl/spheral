text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "DEM/SolidBoundary/CylinderSolidBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CylinderSolidBoundary< Dim< %(ndim)s > >;
}
"""
