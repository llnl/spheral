//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Boundary/findNodesTouchingThroughPlanes.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template
std::vector<size_t> findNodesTouchingThroughPlanes(const NodeList<Dim<1>>& nodeList,
                                                   const GeomPlane<Dim<1>>& enterPlane,
                                                   const GeomPlane<Dim<1>>& exitPlane,
                                                   const double hmultiplier);
#endif

#if defined(SPHERAL_ENABLE_2D)
template
std::vector<size_t> findNodesTouchingThroughPlanes(const NodeList<Dim<2>>& nodeList,
                                                   const GeomPlane<Dim<2>>& enterPlane,
                                                   const GeomPlane<Dim<2>>& exitPlane,
                                                   const double hmultiplier);
#endif

#if defined(SPHERAL_ENABLE_3D)
template
std::vector<size_t> findNodesTouchingThroughPlanes(const NodeList<Dim<3>>& nodeList,
                                                   const GeomPlane<Dim<3>>& enterPlane,
                                                   const GeomPlane<Dim<3>>& exitPlane,
                                                   const double hmultiplier);
#endif

}
