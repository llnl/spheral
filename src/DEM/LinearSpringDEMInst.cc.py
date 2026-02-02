text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "DEM/LinearSpringDEM.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class LinearSpringDEM< Dim< %(ndim)s > >;
}
"""
