text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Neighbor/NestedGridNeighbor.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

template class NestedGridNeighbor< Dim< %(ndim)s > >;

}
"""
