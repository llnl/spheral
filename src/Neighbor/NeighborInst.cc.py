text = """
#include "Neighbor/Neighbor.cc"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class Neighbor<Dim<%(ndim)s>>;
}
"""
