text = """
#include "Neighbor/PythonNeighbor.cc"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class PythonNeighbor<Dim<%(ndim)s>>;
  template class PythonNeighborCallView<Dim<%(ndim)s>>;
}
"""
