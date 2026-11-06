text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Boundary/findNodesTouchingThroughPlanes.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template
  std::vector<size_t> 
  findNodesTouchingThroughPlanes(const NodeList<Dim<%(ndim)s>>&,
                                 const GeomPlane<Dim<%(ndim)s>>&,
                                 const GeomPlane<Dim<%(ndim)s>>&,
                                 const double,
                                 const double);
                         
}
"""
