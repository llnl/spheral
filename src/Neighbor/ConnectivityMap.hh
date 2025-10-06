//---------------------------------Spheral++----------------------------------//
// ConnectivityMap
//
// Stores the full set of significant neighbors for a set of NodeLists.
//
// Created by J. Michael Owen, Sun Oct 30 15:36:33 PST 2005
//----------------------------------------------------------------------------//
#ifndef _Spheral_NeighborSpace_ConnectivityMap_hh_
#define _Spheral_NeighborSpace_ConnectivityMap_hh_

#include "Utilities/KeyTraits.hh"
#include "Utilities/NodeCoupling.hh"
#include "Field/FieldList.hh"
#include "Neighbor/NodePairList.hh"

#include <vector>
#include <unordered_map>
#include <map>
#include <memory>

namespace Spheral {

template<typename Dimension> class NodeList;
template<typename Dimension> class FluidNodeList;
template<typename Dimension> class Boundary;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
class ConnectivityMap {
public:
  //--------------------------- Public Interface ---------------------------//
  using const_iterator = std::vector<int>::const_iterator;
  using NodeCouplingPtr = std::shared_ptr<NodeCoupling>;
  using Vector = typename Dimension::Vector;

  // Constructors, destructor.
  ConnectivityMap();
  ~ConnectivityMap();

  template<typename NodeListIterator>
  ConnectivityMap(const NodeListIterator& begin,
                  const NodeListIterator& end,
                  const bool buildGhostConnectivity,
                  const bool buildOverlapConnectivity,
                  const bool buildIntersectionConnectivity);

  // Rebuild for a given set of NodeLists.
  template<typename NodeListIterator>
  void rebuild(const NodeListIterator& begin, 
               const NodeListIterator& end, 
               const bool computeGhostConnectivity,
               const bool buildOverlapConnectivity,
               const bool buildIntersectionConnectivity);

  // Patch the connectivity information:
  // flags   -- (0,1): 0 => node deleted, 1 => node preserved
  // old2new -- maps old -> new node indices.
  void patchConnectivity(const FieldList<Dimension, size_t>& flags,
                         const FieldList<Dimension, size_t>& old2new);

  // Remove connectivity between neighbors.
  // Note this method assumes neighbor info is symmetric, and removes the pair connectivity for each
  // member of a pair (maintaining symmetry).
  void removeConnectivity(const FieldList<Dimension, std::vector<std::vector<int>>>& neighborsToCut);

  // Remove pairs based on a function
  void removeConnectivity(std::function<bool(int, int, int, int)> excludePairs);

  // Are we computing neighbors for ghosts?
  bool buildGhostConnectivity() const;

  // Do we compute overlap connectivity?
  bool buildOverlapConnectivity() const;

  // Do we compute intersection connectivity?
  bool buildIntersectionConnectivity() const;

  // Get the set of NodeLists.
  const std::vector<const NodeList<Dimension>*>& nodeLists() const;

  // The set of pairs
  const NodePairList& nodePairList() const;
  std::shared_ptr<NodePairList> nodePairListPtr() const;

  // A functor to specify the coupling between nodes
  NodeCoupling& coupling();
  const NodeCoupling& coupling() const;
  void coupling(NodeCouplingPtr couplingPtr);

  //............................................................................
  // Get the set of neighbors for the given (internal!) node in the given NodeList.
  const std::vector< std::vector<int> >&
  connectivityForNode(const NodeList<Dimension>* nodeListPtr,
                      const int nodeID) const;

  // Same as above, just referencing the NodeList by an integer index.
  const std::vector< std::vector<int> >&
  connectivityForNode(const int nodeListID,
                      const int nodeID) const;

  //............................................................................
  // Get the pre-computed intersection connectivity for points if it was requested.
  // Note, this is different than what we expect for overlap connectivity: in this
  // method the intersection points (k) are all points that points (i,j) have in
  // common when (i,j) are ALSO neighbors.  Overlap connectivity may exist for
  // (i,j) even if (i,j) are not neighbors, and this set will miss such points.
  const std::vector<std::vector<int>>&
  intersectionConnectivity(const NodePairIdxType& pair) const;

  //............................................................................
  // Note the following two methods return the points we have neighbors in common with,
  // not the common neighbors.  You need to query ConnectivityMap::connectivityIntersectionForNodes
  // to get the overlapping set of points.
  // Get the set of neighbors we have overlap with (common neighbors).
  const std::vector< std::vector<int> >&
  overlapConnectivityForNode(const NodeList<Dimension>* nodeListPtr,
                             const int nodeID) const;

  // Same as above, just referencing the NodeList by an integer index.
  const std::vector< std::vector<int> >&
  overlapConnectivityForNode(const int nodeListID,
                             const int nodeID) const;

  //............................................................................
  // Compute the common neighbors for a pair of nodes.  Note this method 
  // returns by value since this information is not stored by ConnectivityMap.
  std::vector< std::vector<int> >
  connectivityIntersectionForNodes(const int nodeListi, const int i,
                                   const int nodeListj, const int j,
                                   const FieldList<Dimension, Vector>& position = FieldList<Dimension, Vector>()) const;

  // Compute the union of neighbors for a pair of nodes.  Note this method 
  // returns by value since this information is not stored by ConnectivityMap.
  std::vector< std::vector<int> >
  connectivityUnionForNodes(const int nodeListi, const int i,
                            const int nodeListj, const int j) const;

  // Compute the number of neighbors for the given node.
  size_t numNeighborsForNode(const NodeList<Dimension>* nodeListPtr,
                             const int nodeID) const;

  size_t numNeighborsForNode(const int nodeListID,
                             const int nodeID) const;

  // Compute the number of overlap neighbors for the given node.
  size_t numOverlapNeighborsForNode(const NodeList<Dimension>* nodeListPtr,
                                    const int nodeID) const;

  size_t numOverlapNeighborsForNode(const int nodeListID,
                                    const int nodeID) const;
  
  // Return the connectivity in terms of global node IDs.
  std::map<int, std::vector<int> > globalConnectivity(std::vector<Boundary<Dimension>*>& boundaries) const;

  // Function to determine if given node information (i and j), if the 
  // pair should already have been calculated by iterating over each
  // others neighbors.
  bool calculatePairInteraction(const int nodeListi, const int i, 
                                const int nodeListj, const int j,
                                const int firstGhostNodej) const;

  // Provide iterator interface for walking the nodes in a NodeList
  // in order to maintain domain decomposition independence when 
  // desired.
  const_iterator begin(const int nodeList) const;
  const_iterator end(const int nodeList) const;

  // Return the number of nodes we should walk for the given NodeList.
  size_t numNodes(const int nodeList) const;

  // The ith node (ordered) in the given NodeList.
  int ithNode(const int nodeList, const int index) const;

  // Get the ith NodeList or FluidNodeList.
  const NodeList<Dimension>& nodeList(const int index) const;

  // Return which NodeList index in order the given one would be in our connectivity.
  unsigned nodeListIndex(const NodeList<Dimension>* nodeListPtr) const;

  // Exclude (nodeListi, nodei, nodeListj, nodej) pairs if this returns false
  void setNodePairExclusion(std::function<bool(int, int, int, int)> excludePairs);

  // Check that the internal data structure is valid.
  bool valid() const;
private:
  //--------------------------- Private Interface ---------------------------//
  // The set of NodeLists.
  std::vector<const NodeList<Dimension>*> mNodeLists;

  // Are we building ghost and/or overlap connectivity?
  bool mBuildGhostConnectivity, mBuildOverlapConnectivity, mBuildIntersectionConnectivity;

  // The full connectivity map.  This might be quite large!
  // [offset[NodeList] + nodeID] [NodeListID] [neighborIndex]
  using ConnectivityStorageType = std::vector<std::vector<std::vector<int>>>;
  std::vector<int> mOffsets;
  ConnectivityStorageType mConnectivity;

  // List of Node conncetion pairs.
  std::shared_ptr<NodePairList> mNodePairListPtr;

  // Same for overlap connectivity.
  ConnectivityStorageType mOverlapConnectivity;

  // The set of node indices per Nodelist in order for traversal.
  std::vector<std::vector<int>> mNodeTraversalIndices;

  // The set of keys we may compute for each node.
  using Key = typename KeyTraits::Key;
  FieldList<Dimension, Key> mKeys;

  // The coupling functor for coupling between points.
  NodeCouplingPtr mCouplingPtr;

  // The connectivity intersections, if we're keeping them.
  using IntersectionConnectivityContainer = std::unordered_map<NodePairIdxType, std::vector<std::vector<int>>>;
  IntersectionConnectivityContainer mIntersectionConnectivity;

  // Exclude pairs matching this function
  std::function<bool(int, int, int, int)> mExcludePairs;
   
  // Internal method to fill in the connectivity, once the set of NodeLists 
  // is determined.
  void computeConnectivity();

  // No default constructor, copying, or assignment.
  ConnectivityMap(const ConnectivityMap&);
  ConnectivityMap& operator=(const ConnectivityMap&);
};

}

#include "ConnectivityMapInline.hh"

#endif

