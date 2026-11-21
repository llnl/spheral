text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/NestedGridRedistributeNodes.cc"

namespace Spheral {
  template class NestedGridRedistributeNodes< Dim< %(ndim)s > >;
}

"""
