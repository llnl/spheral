text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "DEM/DEMBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class DEMBase< Dim< %(ndim)s > >;
}
"""
