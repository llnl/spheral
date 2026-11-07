text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/NestedGridDistributedBoundary.cc"

namespace Spheral {
  template class NestedGridDistributedBoundary< Dim< %(ndim)s > >;
}
"""
