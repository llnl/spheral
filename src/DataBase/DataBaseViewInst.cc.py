text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase/DataBaseView.cc"

namespace Spheral {
  template class DataBaseView<Dim<%(ndim)s>>;
}
"""
