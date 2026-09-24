//---------------------------------Spheral++----------------------------------//
// Neighbor -- Abstract interface base class for the Neighbor objects.
//
// Created by J. Michael Owen, Sun Nov 12 10:33:55 2000
//----------------------------------------------------------------------------//
#include "Neighbor/Neighbor.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"
#include "Utilities/DBC.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/NodeListRegistrar.hh"
#include "Utilities/testBoxIntersection.hh"

#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given NodeList and search type.
//------------------------------------------------------------------------------
template<typename Dimension>
Neighbor<Dimension>::
Neighbor(NodeList<Dimension>& nodeList,
         const NeighborSearchType searchType,
         const Scalar kernelExtent):
  NeighborBase<Dimension>(searchType, kernelExtent),
  mNodeListPtr(&nodeList),
  mNodeExtent("Node Extent", nodeList) {
  mNodeListPtr->registerNeighbor(*this);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
Neighbor<Dimension>::
~Neighbor() {
  if (mNodeListPtr != nullptr) mNodeListPtr->unregisterNeighbor();
}

//------------------------------------------------------------------------------
// Set the node list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Neighbor<Dimension>::
nodeList(NodeList<Dimension>& nodeList) {
  mNodeListPtr = &nodeList;
  mNodeExtent.setNodeList(nodeList);
}

//------------------------------------------------------------------------------
// Unregister the node list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Neighbor<Dimension>::
unregisterNodeList() {
  mNodeListPtr = nullptr;
}

}
