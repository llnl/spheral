text = """
#include "Neighbor/NeighborView.cc"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class NeighborView<Dim<%(ndim)s>>;
}
"""
