//---------------------------------Spheral++----------------------------------//
// Neighbor -- Abstract interface base class for the Neighbor objects.
//
// Created by J. Michael Owen, Sun Nov 12 10:33:55 2000
//----------------------------------------------------------------------------//
#include "Neighbor.hh"

#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"
#include "Utilities/DBC.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/NodeListRegistrar.hh"
#include "Utilities/testBoxIntersection.hh"

#include <algorithm>
using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given NodeList and search type.
//------------------------------------------------------------------------------
template<typename Dimension>
Neighbor<Dimension>::
Neighbor(NodeList<Dimension>& nodeList,
         const NeighborSearchType searchType,
         const double kernelExtent):
  mSearchType(searchType),
  mKernelExtent(kernelExtent),
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
// Get the type of search.
//------------------------------------------------------------------------------
template<typename Dimension>
NeighborSearchType
Neighbor<Dimension>::
neighborSearchType() const {
  return mSearchType;
}

//------------------------------------------------------------------------------
// Set the type of search.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Neighbor<Dimension>::
neighborSearchType(NeighborSearchType searchType) {
  mSearchType = searchType;
}

//------------------------------------------------------------------------------
// Return the field of node extents.
//------------------------------------------------------------------------------
template<typename Dimension>
const Field<Dimension, typename Dimension::Vector>&
Neighbor<Dimension>::
nodeExtentField() const {
  return mNodeExtent;
}

// Allow read/write access to the node extent Field for descendent classes.
template<typename Dimension>
Field<Dimension, typename Dimension::Vector>&
Neighbor<Dimension>::
accessNodeExtentField() {
  return mNodeExtent;
}

//------------------------------------------------------------------------------
// Access the node list.
//------------------------------------------------------------------------------
template<typename Dimension>
const NodeList<Dimension>&
Neighbor<Dimension>::
nodeList() const {
  CHECK(mNodeListPtr != nullptr);
  return *mNodeListPtr;
}

template<typename Dimension>
const NodeList<Dimension>*
Neighbor<Dimension>::
nodeListPtr() const {
  return mNodeListPtr;
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

template<typename Dimension>
void
Neighbor<Dimension>::
nodeListPtr(NodeList<Dimension>* nodeListPtr) {
  CHECK(nodeListPtr);
  mNodeListPtr = nodeListPtr;
  mNodeExtent.setNodeList(*nodeListPtr);
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

//------------------------------------------------------------------------------
// Determine the extent of the given node's smoothing tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Vector
Neighbor<Dimension>::
nodeExtent(int nodeID) const {
  CHECK(nodeID >= 0 and nodeID < (int)nodeList().numNodes());
  return HExtent(nodeList().Hfield()(nodeID), kernelExtent());
}

//------------------------------------------------------------------------------
// Force the node extent field to be computed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Neighbor<Dimension>::
setNodeExtents() {
  for (auto nodeID = 0u; nodeID < nodeList().numNodes(); ++nodeID) {
    mNodeExtent(nodeID) = nodeExtent(nodeID);
  }
}

template<typename Dimension>
void
Neighbor<Dimension>::
setNodeExtents(const vector<int>& nodeIDs) {
  for (typename vector<int>::const_iterator nodeIDItr = nodeIDs.begin();
       nodeIDItr < nodeIDs.end();
       ++nodeIDItr) {
    CHECK(*nodeIDItr >= 0 and *nodeIDItr < (int)nodeList().numNodes());
    mNodeExtent(*nodeIDItr) = nodeExtent(*nodeIDItr);
  }
}

template<typename Dimension>
void
Neighbor<Dimension>::
setInternalNodeExtents() {
  for (auto nodeID = 0u; nodeID < nodeList().numInternalNodes(); ++nodeID) {
    mNodeExtent(nodeID) = nodeExtent(nodeID);
  }
}

template<typename Dimension>
void
Neighbor<Dimension>::
setGhostNodeExtents() {
  for (int nodeID = nodeList().firstGhostNode(); nodeID < (int)nodeList().numNodes(); ++nodeID) {
    mNodeExtent(nodeID) = nodeExtent(nodeID);
  }
}

//------------------------------------------------------------------------------
// Set the master list of nodes based on a particular node index.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Neighbor<Dimension>::
setMasterList(int nodeID,
              std::vector<int>& masterList,
              std::vector<int>& coarseNeighbors,
              const bool ghostConnectivity) const {
  CHECK(valid());
  CHECK(nodeID >= 0 and nodeID < (int)nodeList().numInternalNodes());
  const auto& positions = nodeList().positions();
  const auto& Hfield = nodeList().Hfield();
  this->setMasterList(positions(nodeID), Hfield(nodeID), masterList, coarseNeighbors, ghostConnectivity);
}

//------------------------------------------------------------------------------
// Set the refined list of potential neighbors based on a particular node index.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Neighbor<Dimension>::
setRefineNeighborList(int nodeID,
                      const std::vector<int>& coarseNeighbors,
                      std::vector<int>& refineNeighbors) const {
  CHECK(valid());
  CHECK(nodeID >= 0 and nodeID < (int)nodeList().numInternalNodes());
  const auto& positions = nodeList().positions();
  const auto& Hfield = nodeList().Hfield();
  this->setRefineNeighborList(positions(nodeID), Hfield(nodeID), coarseNeighbors, refineNeighbors);
}

//------------------------------------------------------------------------------
// General precull routine.  Select nodes based on a range of positions and
// extents.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<int>
Neighbor<Dimension>::
precullList(const Vector& minMasterPosition, const Vector& maxMasterPosition,
            const Vector& minMasterExtent, const Vector& maxMasterExtent,
            const vector<int>& coarseList) const {
  const auto n = coarseList.size();

  // Empty vector to accumulate the result in.
  vector<int> cullList;

  // Get a reference to the node positions and node extent field.
  const auto& positions = nodeList().positions();
  const auto& nodeExtents = nodeExtentField();

  // What kind of preculling we're doing determines the applied test.
  if (neighborSearchType() == NeighborSearchType::GatherScatter) {
    
    // Gather-Scatter.
// #pragma omp parallel
    // {
    //   vector<int> cullList_private;
// #pragma omp for nowait
      for (auto k = 0u; k < n; ++k) {
        const auto  j = coarseList[k];
        const auto& nodePosition = positions(j);
        const auto  minNodeExtent = nodePosition - nodeExtents(j);
        const auto  maxNodeExtent = nodePosition + nodeExtents(j);
        if (testPointInBox(nodePosition,
                           minMasterExtent,
                           maxMasterExtent) or
            testBoxIntersection(minMasterPosition,
                                maxMasterPosition,
                                minNodeExtent,
                                maxNodeExtent)) {
          cullList.push_back(j);
        }
      }

// #pragma omp critical
//       std::copy(cullList_private.begin(), cullList_private.end(), std::back_inserter(cullList));
//     }

  } else if (neighborSearchType() == NeighborSearchType::Gather) {

    // Gather.
    for (typename vector<int>::const_iterator coarseItr = coarseList.begin();
         coarseItr != coarseList.end();
         ++coarseItr) {
      const auto  j = *coarseItr;
      const auto& nodePosition = positions(j);
      const auto  gatherTest = testPointInBox(nodePosition,
                                              minMasterExtent,
                                              maxMasterExtent);
      if (gatherTest) cullList.push_back(j);
    }

  } else {

    // Scatter.
    CHECK(neighborSearchType() == NeighborSearchType::Scatter);
    for (typename vector<int>::const_iterator coarseItr = coarseList.begin();
         coarseItr != coarseList.end();
         ++coarseItr) {
      const auto  j = *coarseItr;
      const auto& nodePosition = positions(j);
      const auto  minNodeExtent = nodePosition - nodeExtents(j);
      const auto  maxNodeExtent = nodePosition + nodeExtents(j);
      const auto  scatterTest = testBoxIntersection(minMasterPosition,
                                                    maxMasterPosition,
                                                    minNodeExtent,
                                                    maxNodeExtent);
      if (scatterTest) cullList.push_back(j);
    }

  }

  ENSURE(cullList.size() <= coarseList.size());
  return cullList;
}

// //------------------------------------------------------------------------------
// // Cull the local (to this NodeList) neighbor info based on the current master
// // state.
// // *NOTE* -- this is not safe to do when you want to use this neighbor info 
// // with different NodeLists!
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void
// Neighbor<Dimension>::
// precullForLocalNodeList() {

//   // Grab the state.
//   const Field<Dimension, Vector>& r = this->nodeList().positions();
//   const Field<Dimension, Vector>& extent = nodeExtentField();

//   // Find the min/max master node positions and extents.
//   Vector minMasterPosition = DBL_MAX;
//   Vector maxMasterPosition = -DBL_MAX;
//   Vector minMasterExtent = DBL_MAX;
//   Vector maxMasterExtent = -DBL_MAX;
//   for (const_iterator masterItr = masterBegin();
//        masterItr != masterEnd();
//        ++masterItr) {
//     const int i = *masterItr;
//     const Vector& ri = r(i);
//     const Vector minExtent = ri - extent(i);
//     const Vector maxExtent = ri + extent(i);
//     for (int j = 0; j != Dimension::nDim; ++j) {
//       minMasterPosition(j) = min(minMasterPosition(j), ri(j));
//       maxMasterPosition(j) = max(maxMasterPosition(j), ri(j));
//       minMasterExtent(j) = min(minMasterExtent(j), minExtent(j));
//       maxMasterExtent(j) = max(maxMasterExtent(j), maxExtent(j));
//     }
//   }

//   // Now use this info to precull the coarse neighbor list.
//   *mCoarseNeighborListPtr = precullList(minMasterPosition, maxMasterPosition,
//                                         minMasterExtent, maxMasterExtent,
//                                         *mCoarseNeighborListPtr);
// }

//------------------------------------------------------------------------------
// Provide a basic test to determine whether the neighbor base is in a valid
// (i.e., ready to use) state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
Neighbor<Dimension>::
valid() const {
  return (kernelExtent() > 0.0);
}

// //------------------------------------------------------------------------------
// // Set the global bounding box.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void
// Neighbor<Dimension>::
// setBoundingBox() {
//   const NodeListRegistrar<Dimension>& registrar = NodeListRegistrar<Dimension>::instance();
//   FieldList<Dimension, Vector> positions(FieldListBase::Reference);
//   for (typename NodeListRegistrar<Dimension>::const_fluid_iterator itr = registrar.fluidBegin();
//        itr != registrar.fluidEnd();
//        ++itr) {
//     positions.appendField((**itr).positions());
//   }
//   globalBoundingBox(positions, mXmin, mXmax, false);
// }

}
