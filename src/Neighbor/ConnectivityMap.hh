//---------------------------------Spheral++----------------------------------//
// ConnectivityMap
//
// Stores the full set of significant neighbors for a set of NodeLists.
//
// Created by J. Michael Owen, Sun Oct 30 15:36:33 PST 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_ConnectivityMap__
#define __Spheral_ConnectivityMap__

#include "Utilities/KeyTraits.hh"
#include "Utilities/NodeCoupling.hh"
#include "Field/FieldList.hh"
#include "Neighbor/NodePairList.hh"
#include "Neighbor/NeighborVariant.hh"

#include <vector>
#include <unordered_map>
#include <map>
#include <memory>
#include <algorithm>

namespace Spheral {

template<typename Dimension> class NodeList;
template<typename Dimension> class FluidNodeList;
template<typename Dimension> class DataBase;
template<typename Dimension> class Boundary;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
class ConnectivityMap {
public:
  //--------------------------- Public Interface ---------------------------//
  using const_iterator = std::vector<size_t>::const_iterator;
  using NodeCouplingPtr = std::shared_ptr<NodeCoupling>;
  using Vector = typename Dimension::Vector;
  using SymTensor = typename Dimension::SymTensor;

  using NodeID = std::pair<size_t, size_t>;  // (NodeList ID, node ID)

  // Constructors, destructor.
  ConnectivityMap(const DataBase<Dimension>& dataBase,
                  const bool buildGhostConnectivity,
                  const bool buildOverlapConnectivity,
                  const bool buildIntersectionConnectivity);
  ~ConnectivityMap() = default;

  // Rebuild for a given set of NodeLists.
  void rebuild(const bool computeGhostConnectivity,
               const bool buildOverlapConnectivity,
               const bool buildIntersectionConnectivity);

  // Patch the connectivity information:
  // flags   -- (0,1): 0 => node deleted, 1 => node preserved
  // old2new -- maps old -> new node indices.
  void patchConnectivity(const FieldList<Dimension, size_t>& flags,
                         const FieldList<Dimension, size_t>& old2new);

  // Access our state
  const DataBase<Dimension>& dataBase()                      const { return mDataBase; }

  // Are we computing neighbors for ghosts?
  bool buildGhostConnectivity()                              const { return mBuildGhostConnectivity; }

  // Do we compute overlap connectivity?
  bool buildOverlapConnectivity()                            const { return mBuildOverlapConnectivity; }

  // Do we compute intersection connectivity?
  bool buildIntersectionConnectivity()                       const { return mBuildIntersectionConnectivity; }

  // Get the set of NodeLists.
  const std::vector<NodeList<Dimension>*>& nodeLists()       const { return mDataBase.nodeListPtrs(); }

  // The set of pairs
  const NodePairList& nodePairList()                         const { REQUIRE(mNodePairListPtr); return *mNodePairListPtr; }
  std::shared_ptr<NodePairList> nodePairListPtr()            const { return mNodePairListPtr; }

  // A functor to specify the coupling between nodes
  NodeCoupling& coupling()                                         { return *mCouplingPtr; }
  const NodeCoupling& coupling()                             const { return *mCouplingPtr; }
  void coupling(NodeCouplingPtr couplingPtr)                       { mCouplingPtr = couplingPtr; }

  //............................................................................
  // Get the set of neighbors for the given (internal!) node in the given NodeList.
  const std::vector<std::vector<size_t>>&
  connectivityForNode(const NodeList<Dimension>* nodeListPtr,
                      const size_t nodeID) const;

  // Same as above, just referencing the NodeList by an integer index.
  const std::vector<std::vector<size_t>>&
  connectivityForNode(const size_t nodeListID,
                      const size_t nodeID) const;

  //............................................................................
  // Get the pre-computed intersection connectivity for points if it was requested.
  // Note, this is different than what we expect for overlap connectivity: in this
  // method the intersection points (k) are all points that points (i,j) have in
  // common when (i,j) are ALSO neighbors.  Overlap connectivity may exist for
  // (i,j) even if (i,j) are not neighbors, and this set will miss such points.
  const std::vector<std::vector<size_t>>&
  intersectionConnectivity(const NodePairIdxType& pair) const;

  //............................................................................
  // Note the following two methods return the points we have neighbors in common with,
  // not the common neighbors.  You need to query ConnectivityMap::connectivityIntersectionForNodes
  // to get the overlapping set of points.
  // Get the set of neighbors we have overlap with (common neighbors).
  const std::vector<std::vector<size_t>>&
  overlapConnectivityForNode(const NodeList<Dimension>* nodeListPtr,
                             const size_t nodeID) const;

  // Same as above, just referencing the NodeList by an integer index.
  const std::vector<std::vector<size_t>>&
  overlapConnectivityForNode(const size_t nodeListID,
                             const size_t nodeID) const;

  //............................................................................
  // Compute the common neighbors for a pair of nodes.  Note this method 
  // returns by value since this information is not stored by ConnectivityMap.
  std::vector<std::vector<size_t>>
  connectivityIntersectionForNodes(const size_t nodeListi, const size_t i,
                                   const size_t nodeListj, const size_t j,
                                   const FieldList<Dimension, Vector>& position = FieldList<Dimension, Vector>()) const;

  // Compute the union of neighbors for a pair of nodes.  Note this method 
  // returns by value since this information is not stored by ConnectivityMap.
  std::vector<std::vector<size_t>>
  connectivityUnionForNodes(const size_t nodeListi, const size_t i,
                            const size_t nodeListj, const size_t j) const;

  // Compute the number of neighbors for the given node.
  size_t numNeighborsForNode(const NodeList<Dimension>* nodeListPtr,
                             const size_t nodeID) const;

  size_t numNeighborsForNode(const size_t nodeListID,
                             const size_t nodeID) const;

  // Compute the number of overlap neighbors for the given node.
  size_t numOverlapNeighborsForNode(const NodeList<Dimension>* nodeListPtr,
                                    const size_t nodeID) const;

  size_t numOverlapNeighborsForNode(const size_t nodeListID,
                                    const size_t nodeID) const;
  
  // Return the connectivity in terms of global node IDs.
  std::map<size_t, std::vector<size_t> > globalConnectivity(std::vector<Boundary<Dimension>*>& boundaries) const;

  // Function to determine if given node information (i and j), if the 
  // pair should already have been calculated by iterating over each
  // others neighbors.
  bool calculatePairInteraction(const size_t nodeListi, const size_t i, 
                                const size_t nodeListj, const size_t j,
                                const size_t firstGhostNodej) const;

  // Provide iterator interface for walking the nodes in a NodeList
  // in order to maintain domain decomposition independence when 
  // desired.
  const_iterator begin(const size_t nodeList) const;
  const_iterator end(const size_t nodeList) const;

  // Return the number of nodes we should walk for the given NodeList.
  size_t numNodes(const size_t nodeList) const;

  // The ith node (ordered) in the given NodeList.
  size_t ithNode(const size_t nodeList, const size_t index) const;

  // Get the ith NodeList or FluidNodeList.
  const NodeList<Dimension>& nodeList(const size_t index) const;

  // Check that the internal data structure is valid.
  bool valid() const;

  // Forbidden methods
  ConnectivityMap() = delete;
  ConnectivityMap(const ConnectivityMap&) = delete;
  ConnectivityMap& operator=(const ConnectivityMap&) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  const DataBase<Dimension>& mDataBase;

  // Are we building ghost and/or overlap connectivity?
  bool mBuildGhostConnectivity, mBuildOverlapConnectivity, mBuildIntersectionConnectivity;

  // The full connectivity map.  This might be quite large!
  // [offset[NodeList] + nodeID] [NodeListID] [neighborIndex]
  // Note this is only built as needed, we assume only the NodePairList by default
  using ConnectivityStorageType = std::vector<std::vector<std::vector<size_t>>>;
  std::vector<size_t> mOffsets;
  ConnectivityStorageType mConnectivity;

  // List of Node connection pairs.
  std::shared_ptr<NodePairList> mNodePairListPtr;

  // Same for overlap connectivity.
  ConnectivityStorageType mOverlapConnectivity;

  // The set of node indices per Nodelist in order for traversal.
  std::vector<std::vector<size_t>> mNodeTraversalIndices;

  // The set of keys we may compute for each node.
  using Key = typename KeyTraits::Key;
  FieldList<Dimension, Key> mKeys;

  // The number of neighbors per node
  FieldList<Dimension, size_t> mNumNeighbors;

  // The coupling functor for coupling between points.
  NodeCouplingPtr mCouplingPtr;

  // The connectivity intersections, if we're keeping them.
  using IntersectionConnectivityContainer = std::unordered_map<NodePairIdxType, std::vector<std::vector<size_t>>>;
  IntersectionConnectivityContainer mIntersectionConnectivity;

  // Internal method to fill in the connectivity, once the set of NodeLists 
  // is determined.
  void computeConnectivity();
  template<typename NeighborType> void computeConnectivityImpl(const NeighborType& neighbor);
  void buildPerPointConnectivity();
};

}

#include "ConnectivityMapInline.hh"

#endif

