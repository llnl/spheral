//---------------------------------Spheral++----------------------------------//
// Find the set of nodes that see through a pair of planes.
//
// Created by JMO, Wed Oct 16 10:00:14 PDT 2019
//
// Modified by:
//----------------------------------------------------------------------------//
#ifndef __Spheral_findNodesTouchingThroughPlanes__
#define __Spheral_findNodesTouchingThroughPlanes__

#include "Geometry/GeomPlane.hh"
#include "NodeList/NodeList.hh"

#include <vector>

namespace Spheral {
  
template<typename Dimension>
std::vector<size_t>
findNodesTouchingThroughPlanes(const NodeList<Dimension>& nodeList,
                               const GeomPlane<Dimension>& enterPlane,
                               const GeomPlane<Dimension>& exitPlane,
                               const double hmultiplier = 1.0);

}

#endif
