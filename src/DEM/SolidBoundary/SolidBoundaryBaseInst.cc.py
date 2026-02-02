text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "DEM/SolidBoundary/SolidBoundaryBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SolidBoundaryBase< Dim< %(ndim)s > >;
}
"""
