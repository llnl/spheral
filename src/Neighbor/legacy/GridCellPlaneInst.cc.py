text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Neighbor/GridCellPlane.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class GridCellPlane< Dim< %(ndim)s > >;
}
"""
