//---------------------------------Spheral++----------------------------------//
// Boundary -- Abstract base class for the boundary condition classes.
//
// Created by JMO, Wed Feb 16 21:01:02 PST 2000
//----------------------------------------------------------------------------//

#include <algorithm>

#include "Boundary.hh"
#include "DataBase/DataBase.hh"
#include "NodeList/NodeList.hh"
#include "Field/FieldList.hh"
#include "Utilities/removeElements.hh"
#include "RK/RKCoefficients.hh"

#include "Utilities/DBC.hh"

using std::vector;
using std::map;

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
Boundary<Dimension>::Boundary():
  mBoundaryNodes() {
}

//------------------------------------------------------------------------------
// Default method for handing ghost nodes for Fields -- copy control->ghost values
//------------------------------------------------------------------------------
template<typename Dimension>
void
Boundary<Dimension>::applyGhostBoundary(FieldBase<Dimension>& fieldBase) const {
  const auto& boundaryNodes = accessBoundaryNodes(const_cast<NodeList<Dimension>&>(fieldBase.nodeList()));
  CHECK(boundaryNodes.ghostNodes.size() == boundaryNodes.controlNodes.size());
  if (not boundaryNodes.ghostNodes.empty()) fieldBase.copyElements(boundaryNodes.controlNodes,
                                                                   boundaryNodes.ghostNodes);
}

//------------------------------------------------------------------------------
// Default method to set ghost nodes for all NodeLists in DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Boundary<Dimension>::setAllGhostNodes(DataBase<Dimension>& dataBase) {
  for (auto nodeListItr = dataBase.nodeListBegin();
       nodeListItr < dataBase.nodeListEnd();
       ++nodeListItr) setGhostNodes(**nodeListItr);
}

//------------------------------------------------------------------------------
// Default method to set violation nodes for all NodeLists in DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Boundary<Dimension>::setAllViolationNodes(DataBase<Dimension>& dataBase) {
  for (typename DataBase<Dimension>::NodeListIterator 
	 nodeListItr = dataBase.nodeListBegin();
       nodeListItr < dataBase.nodeListEnd();
       ++nodeListItr) {
    setViolationNodes(**nodeListItr);
  }

  BEGIN_CONTRACT_SCOPE
  // Make sure the boundary knows about each nodelist in the database.
  for (typename DataBase<Dimension>::ConstNodeListIterator 
       i = dataBase.nodeListBegin(); i != dataBase.nodeListEnd(); ++i)
  {
    ENSURE(mBoundaryNodes.find(*i) != mBoundaryNodes.end());
  } // end for
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Cull out inactive ghost nodes based on a FieldList of flags.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Boundary<Dimension>::cullGhostNodes(const FieldList<Dimension, size_t>& flagSet,
                                    FieldList<Dimension, size_t>& old2newIndexMap,
                                    vector<size_t>& numNodesRemoved) {
  const auto& nodeListPtrs = flagSet.nodeListPtrs();
  REQUIRE(numNodesRemoved.size() == nodeListPtrs.size());
  // auto& registrar = NodeListRegistrar<Dimension>::instance();
  // REQUIRE(numNodesRemoved.size() == registrar.numNodeLists());

  // Walk the NodeLists.
  auto nodeListi = 0;
  for (auto nodeListPtr : nodeListPtrs) {
    // Does the Boundary have entries for this NodeList?
    if (haveNodeList(*nodeListPtr)) {
      auto& boundaryNodes = accessBoundaryNodes(const_cast<NodeList<Dimension>&>(*nodeListPtr));
      CHECK(boundaryNodes.ghostNodes.size() == boundaryNodes.controlNodes.size());
      if (boundaryNodes.ghostNodes.size() > 0) {
        const auto myOldFirstGhostNode = boundaryNodes.ghostNodes[0];
        const auto myNewFirstGhostNode = myOldFirstGhostNode - numNodesRemoved[nodeListi];

        // Grab the flags for this NodeList.
        CHECK(flagSet.haveNodeList(*nodeListPtr));
        const auto& flags = *(flagSet[nodeListi]);

        // Patch up the ghost and control node indices.
        vector<size_t> newGhostNodes, newControlNodes;
        auto newGhostIndex = myNewFirstGhostNode;
        for (size_t k = 0; k < boundaryNodes.ghostNodes.size(); ++k) {
          if (flags(boundaryNodes.ghostNodes[k]) == 1u) {
            newGhostNodes.push_back(newGhostIndex);
            old2newIndexMap(nodeListi, boundaryNodes.ghostNodes[k]) = newGhostIndex;
            ++newGhostIndex;
            // CHECK(boundaryNodes.controlNodes[k] < myOldFirstGhostNode);  // <-- Not true for the ConstantBoundary
            newControlNodes.push_back(old2newIndexMap(nodeListi, boundaryNodes.controlNodes[k]));
            // CHECK(newControlNodes.back() < myNewFirstGhostNode);  // <-- Not true for the ConstantBoundary
          } else {
            ++numNodesRemoved[nodeListi];
          }
        }

        // Update the ghost & control nodes, and the result indicating how many of our ghost nodes 
        // were removed.
        boundaryNodes.ghostNodes = newGhostNodes;
        boundaryNodes.controlNodes = newControlNodes;
        CHECK(boundaryNodes.ghostNodes.size() == boundaryNodes.controlNodes.size());
      }
    }
    ++nodeListi;
  }
}

//------------------------------------------------------------------------------
// Clear out any NodeList information that is currently present.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Boundary<Dimension>::reset(const DataBase<Dimension>&) {
  // Clear our own internal data.
  mBoundaryNodes.clear();
}

//------------------------------------------------------------------------------
// Report the number of ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
int
Boundary<Dimension>::numGhostNodes() const {
  int result = 0;
  for (auto itr = mBoundaryNodes.begin(); itr != mBoundaryNodes.end(); ++itr) {
    result += itr->second.ghostNodes.size();
  }
  return result;
}

//------------------------------------------------------------------------------
// Allow read access to the map of NodeList->BoundaryNodes.
//------------------------------------------------------------------------------
template<typename Dimension>
const map<NodeList<Dimension>*, typename Boundary<Dimension>::BoundaryNodes>&
Boundary<Dimension>::boundaryNodeMap() const {
  return mBoundaryNodes;
}

//------------------------------------------------------------------------------
// Return the list of control nodes for the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
const vector<size_t>&
Boundary<Dimension>::controlNodes(const NodeList<Dimension>& nodeList) const {
  auto itr = mBoundaryNodes.find(const_cast<NodeList<Dimension>*>(&nodeList));
  if (itr == mBoundaryNodes.end()) {
    VERIFY2(false, "Boundary::controlNodes: no entry for NodeList: " << nodeList.name());
  }
  return itr->second.controlNodes;
}

//------------------------------------------------------------------------------
// Return the list of ghost nodes for the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
const vector<size_t>&
Boundary<Dimension>::ghostNodes(const NodeList<Dimension>& nodeList) const {
  auto itr = mBoundaryNodes.find(const_cast<NodeList<Dimension>*>(&nodeList));
  if (itr == mBoundaryNodes.end()) {
    std::cerr << "Number of known NodeLists: " << mBoundaryNodes.size() << std::endl;
    VERIFY2(false, "Boundary::ghostNodes: no entry for NodeList: " << nodeList.name());
  }
  return itr->second.ghostNodes;
}

//------------------------------------------------------------------------------
// Return the list of violation nodes for the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
const vector<size_t>&
Boundary<Dimension>::violationNodes(const NodeList<Dimension>& nodeList) const {
  auto itr = mBoundaryNodes.find(const_cast<NodeList<Dimension>*>(&nodeList));
  if (itr == mBoundaryNodes.end()) {
    VERIFY2(false, "Boundary::violationNodes: no entry for NodeList: " << nodeList.name());
  }
  return itr->second.violationNodes;
}

//------------------------------------------------------------------------------
// Return begin/end iterators for the control nodes of the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
typename vector<size_t>::const_iterator
Boundary<Dimension>::controlBegin(const NodeList<Dimension>& nodeList) const {
  return controlNodes(nodeList).begin();
}

template<typename Dimension>
typename vector<size_t>::const_iterator
Boundary<Dimension>::controlEnd(const NodeList<Dimension>& nodeList) const {
  return controlNodes(nodeList).end();
}

//------------------------------------------------------------------------------
// Return begin/end iterators for the ghost nodes of the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
typename vector<size_t>::const_iterator
Boundary<Dimension>::ghostBegin(const NodeList<Dimension>& nodeList) const {
  return ghostNodes(nodeList).begin();
}

template<typename Dimension>
typename vector<size_t>::const_iterator
Boundary<Dimension>::ghostEnd(const NodeList<Dimension>& nodeList) const {
  return ghostNodes(nodeList).end();
}

//------------------------------------------------------------------------------
// Return begin/end iterators for the violation nodes of the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
typename vector<size_t>::const_iterator
Boundary<Dimension>::violationBegin(const NodeList<Dimension>& nodeList) const {
  return violationNodes(nodeList).begin();
}

template<typename Dimension>
typename vector<size_t>::const_iterator
Boundary<Dimension>::violationEnd(const NodeList<Dimension>& nodeList) const {
  return violationNodes(nodeList).end();
}

//------------------------------------------------------------------------------
// accessBoundaryNodes
//------------------------------------------------------------------------------
template<typename Dimension>
map<NodeList<Dimension>*, typename Boundary<Dimension>::BoundaryNodes>&
Boundary<Dimension>::accessBoundaryNodes() {
  return mBoundaryNodes;
}

template<typename Dimension>
const map<NodeList<Dimension>*, typename Boundary<Dimension>::BoundaryNodes>&
Boundary<Dimension>::accessBoundaryNodes() const {
  return mBoundaryNodes;
}

//------------------------------------------------------------------------------
// accessBoundaryNodes(nodeList)
//------------------------------------------------------------------------------
template<typename Dimension>
typename Boundary<Dimension>::BoundaryNodes&
Boundary<Dimension>::accessBoundaryNodes(NodeList<Dimension>& nodeList) {
  auto itr = mBoundaryNodes.find(&nodeList);
  if (itr == mBoundaryNodes.end()) {
    VERIFY2(false, "Boundary::accessBoundaryNodes: no entry for NodeList: " << nodeList.name());
  }
  return itr->second;
}

template<typename Dimension>
const typename Boundary<Dimension>::BoundaryNodes&
Boundary<Dimension>::accessBoundaryNodes(NodeList<Dimension>& nodeList) const {
  auto itr = mBoundaryNodes.find(&nodeList);
  if (itr == mBoundaryNodes.end()) {
    VERIFY2(false, "Boundary::accessBoundaryNodes: no entry for NodeList: " << nodeList.name());
  }
  return itr->second;
}

//------------------------------------------------------------------------------
// Add the given NodeList to the current set of BoundaryNodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Boundary<Dimension>::addNodeList(NodeList<Dimension>& nodeList) {
  if (mBoundaryNodes.find(&nodeList) == mBoundaryNodes.end()) {
    mBoundaryNodes[&nodeList] = BoundaryNodes();
  }
}

}
