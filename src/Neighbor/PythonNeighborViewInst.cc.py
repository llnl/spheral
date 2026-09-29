text = """
#include "Neighbor/PythonNeighborView.cc"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class PythonNeighborView<Dim<%(ndim)s>>;
}
"""
