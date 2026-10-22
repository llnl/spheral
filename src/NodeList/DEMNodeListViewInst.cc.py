text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/DEMNodeListView.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class DEMNodeListView<Dim<%(ndim)s>>;
}
"""
