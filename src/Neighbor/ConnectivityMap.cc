//---------------------------------Spheral++----------------------------------//
// ConnectivityMap
//
// Stores the full set of significant neighbors for a set of NodeLists.
//
// Created by J. Michael Owen, Sun Oct 30 15:36:33 PST 2005
//----------------------------------------------------------------------------//
#include "ConnectivityMap.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/mortonOrderIndices.hh"
#include "Utilities/PairComparisons.hh"
#include "Utilities/pointDistances.hh"
#include "Utilities/Timer.hh"

#include <algorithm>
#include <ctime>
using std::vector;
using std::map;
using std::string;
using std::pair;

namespace Spheral {

namespace {
//------------------------------------------------------------------------------
// Append v2 to the end of v1
//------------------------------------------------------------------------------
template<typename T>
inline
void
appendSTLvectors(std::vector<T>& v1, std::vector<T>& v2) {
  if (not v2.empty()) {
    v1.reserve(v1.size() + v2.size());
    v1.insert(v1.end(), v2.begin(), v2.end());
  }
}

//------------------------------------------------------------------------------
// How should we compare pairs for sorting?
//------------------------------------------------------------------------------
inline
KeyTraits::Key
hashKeys(const KeyTraits::Key& a, const KeyTraits::Key& b) {
  // Szudzik's function
  // return (a >= b ?
  //         a * a + a + b :
  //         a + b * b);          // where a, b >= 0
  // return a + b;
  return ((KeyTraits::Key(a) << 16) | KeyTraits::Key(b));
}

//------------------------------------------------------------------------------
// Helper to insert into a sorted list of IDs.
//------------------------------------------------------------------------------
template<typename KeyContainer>
inline
bool
insertUnique(const std::vector<size_t>& offsets,
             std::vector<std::vector<std::vector<int>>>& indices,
             const KeyContainer& keys,
             const bool useKeys,
             const int jN1, const int j1,
             const int jN2, const int j2) {
  if (jN1 != jN2 or j1 != j2) {
    auto& overlap = indices[offsets[jN1] + j1][jN2];
    std::vector<int>::iterator itr;
    if (useKeys) {
      itr = std::lower_bound(overlap.begin(), overlap.end(), j2,
                             [&](const int a, const int& b) { return keys(jN2, a) < keys(jN2, b); });
    } else {
      itr = std::lower_bound(overlap.begin(), overlap.end(), j2);
    }
    if (itr == overlap.end() or *itr != j2) {
      overlap.insert(itr, j2);
      return true;
    } else {
      return false;
    }
  }
  return false;
}

//------------------------------------------------------------------------------
// Sort node pairs for domain independent ordering.
//------------------------------------------------------------------------------
template<typename KeyContainer>
inline
void
sortPairs(NodePairList& pairs,
          const KeyContainer& keys) {
  // Start by making sure the pairs themselves are deterministically arranged.
  for (auto& p: pairs) {
    if (keys(p.i_list, p.i_node) > keys(p.j_list, p.j_node)) {
      std::swap(p.i_list, p.j_list);
      std::swap(p.i_node, p.j_node);
    }
  }

  // Now sort the list as a whole.
  std::sort(pairs.begin(), pairs.end(),
            [&](const NodePairIdxType& a, const NodePairIdxType& b) {
              return hashKeys(keys(a.i_list, a.i_node), keys(a.j_list, a.j_node)) <
                hashKeys(keys(b.i_list, b.i_node), keys(b.j_list, b.j_node));
            });
}

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
ConnectivityMap<Dimension>::
ConnectivityMap():
  mNodeLists(),
  mBuildGhostConnectivity(false),
  mBuildOverlapConnectivity(false),
  mBuildIntersectionConnectivity(false),
  mConnectivity(),
  mNodeTraversalIndices(),
  mKeys(FieldStorageType::CopyFields),
  mCouplingPtr(std::make_shared<NodeCoupling>()),
  mIntersectionConnectivity() {
}

//------------------------------------------------------------------------------
// Internal method to build the connectivity for the requested set of NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConnectivityMap<Dimension>::
patchConnectivity(const FieldList<Dimension, size_t>& flags,
                  const FieldList<Dimension, size_t>& old2new) {
  TIME_BEGIN("ConnectivityMap_patch");

  const auto domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();

  // We have to recompute the keys to sort nodes by excluding the 
  // nodes that are being removed.
  const auto numNodeLists = mNodeLists.size();
  if (domainDecompIndependent) {
    for (auto iNodeList = 0u; iNodeList < numNodeLists; ++iNodeList) {
      for (auto i = 0u; i < mNodeLists[iNodeList]->numNodes(); ++i) {
        if (flags(iNodeList, i) == 0) mKeys(iNodeList, i) = KeyTraits::maxKey;
      }
    }
  }

  // Iterate over the traversal ordering and fix it
  for (auto iNodeList = 0u; iNodeList != numNodeLists; ++iNodeList) {
    const auto numNodes = ((domainDecompIndependent or mBuildGhostConnectivity or mBuildOverlapConnectivity) ? 
                           mNodeLists[iNodeList]->numNodes() :
                           mNodeLists[iNodeList]->numInternalNodes());

    vector<size_t> iNodesToKill;
    vector<pair<int, Key>> keys;
#pragma omp parallel
    {
      vector<size_t> iNodesToKill_thread;
      vector<pair<int, Key>> keys_thread;

      // Patch the traversal ordering for this NodeList.
#pragma omp for schedule(dynamic)
      for (auto i = 0u; i < numNodes; ++i) {
        if (flags(iNodeList, i) == 0) {
          iNodesToKill_thread.push_back(i);
        } else {
          if (domainDecompIndependent) keys_thread.push_back(std::make_pair(old2new(iNodeList, i), mKeys(iNodeList, i)));
          mNodeTraversalIndices[iNodeList][i] = old2new(iNodeList, i);
        }
      }

#pragma omp critical
      appendSTLvectors(iNodesToKill, iNodesToKill_thread);
      appendSTLvectors(keys, keys_thread);
    }
    removeElements(mNodeTraversalIndices[iNodeList], iNodesToKill);

    // Recompute the ordering for traversing the nodes.
    {
      const auto numNodes = mNodeTraversalIndices[iNodeList].size();
      if (domainDecompIndependent) {
        std::sort(keys.begin(), keys.end(), ComparePairsBySecondElement<pair<int, Key> >());
#pragma omp parallel for
        for (auto k = 0u; k < numNodes; ++k) {
          mNodeTraversalIndices[iNodeList][k] = keys[k].first;
        }
      } else {
#pragma omp parallel for
        for (auto i = 0u; i < numNodes; ++i) {
          mNodeTraversalIndices[iNodeList][i] = i;
        }
      }
    }
  }
  
  // Patch the node pair structure.
  // Note here we deliberately reallocate the NodePairList, which will invalidate any
  // PairFields pointing at the original pairs.
  REQUIRE(mNodePairListPtr);
  NodePairList& currentPairs = *mNodePairListPtr;
  std::vector<NodePairIdxType> culledPairs;
  culledPairs.reserve(currentPairs.size());
#pragma omp parallel
  {
    std::vector<NodePairIdxType> culledPairs_thread;
    const auto npairs = currentPairs.size();
    culledPairs_thread.reserve(npairs);
#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      const auto iNodeList = currentPairs[k].i_list;
      const auto jNodeList = currentPairs[k].j_list;
      const auto i = currentPairs[k].i_node;
      const auto j = currentPairs[k].j_node;
      if (flags(iNodeList, i) != 0 and flags(jNodeList, j) != 0) {
        culledPairs_thread.push_back(NodePairIdxType(old2new(iNodeList, i), iNodeList,
                                                     old2new(jNodeList, j), jNodeList));
      }
    }
#pragma omp critical
    {
      culledPairs.insert(culledPairs.end(), culledPairs_thread.begin(), culledPairs_thread.end());
    }
  }
  mNodePairListPtr = std::make_shared<NodePairList>(std::move(culledPairs));

  // Sort the NodePairList in order to enforce domain decomposition independence.
  {
    auto& pairs = *mNodePairListPtr;
    if (domainDecompIndependent) {
      sortPairs(pairs, mKeys);
    } else {
      std::sort(pairs.begin(), pairs.end());
    }
  }

  // If needed, rebuild the per-point connectivity
  if (not mConnectivity.empty()) buildPerPointConnectivity();

  // Patch the intersection lists if we're maintaining them
  if (mBuildIntersectionConnectivity) {
    IntersectionConnectivityContainer intersection;
    for (const auto& element: mIntersectionConnectivity) {
      auto pair = element.first;
      const auto& oldintersect = element.second;
      if (flags(pair.i_list, pair.i_node) != 0 and flags(pair.j_list, pair.j_node) != 0) {
        pair.i_node = old2new(pair.i_list, pair.i_node);
        pair.j_node = old2new(pair.j_list, pair.j_node);
        vector<vector<int>> newintersect(numNodeLists);
        for (auto klist = 0u; klist < numNodeLists; ++klist) {
          for (const auto k: oldintersect[klist]) {
            if (flags(klist, k) != 0) newintersect[klist].push_back(old2new(klist, k));
          }
        }
        intersection[pair] = newintersect;
      }
    }
  }

  // You can't check valid yet 'cause the NodeLists have not been resized
  // when we call patch!  The valid method should be checked by whoever called
  // this method after that point.
  //ENSURE(valid());
  TIME_END("ConnectivityMap_patch");
}

//------------------------------------------------------------------------------
// Compute the common neighbors for a pair of nodes.  Note this method 
// returns by value since this information is not stored by ConnectivityMap.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<vector<int> >
ConnectivityMap<Dimension>::
connectivityIntersectionForNodes(const int nodeListi, const int i,
                                 const int nodeListj, const int j,
                                 const FieldList<Dimension, typename Dimension::Vector>& position) const {

  // Pre-conditions.
  TIME_BEGIN("ConnectivityMap_computeIntersectionConnectivity");
  const auto numNodeLists = mNodeLists.size();
  const auto domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
  const auto ghostConnectivity = (mBuildGhostConnectivity or
                                  mBuildOverlapConnectivity or
                                  domainDecompIndependent);
  const auto firstGhostNodei = mNodeLists[nodeListi]->firstGhostNode();
  const auto firstGhostNodej = mNodeLists[nodeListj]->firstGhostNode();
  const bool usePosition = (position.numFields() == numNodeLists);
  REQUIRE(nodeListi < (int)numNodeLists and
          nodeListj < (int)numNodeLists);
  REQUIRE(ghostConnectivity or i < (int)firstGhostNodei or j < (int)firstGhostNodej);
  REQUIRE(position.numFields() == numNodeLists or position.numFields() == 0);

  // Prepare the result.
  vector<vector<int>> result(numNodeLists);

  // If both nodes are internal, we simply intersect their neighbor lists.
  if (ghostConnectivity or (i < (int)firstGhostNodei and j < (int)firstGhostNodej)) {
    const auto& neighborsi = this->connectivityForNode(nodeListi, i);
    const auto& neighborsj = this->connectivityForNode(nodeListj, j);
    CHECK(neighborsi.size() == numNodeLists);
    CHECK(neighborsj.size() == numNodeLists);
    vector<int> neighborsijk;
    Vector posi, posj, b;
    if (usePosition) {
      posi = position(nodeListi, i);
      posj = position(nodeListj, j);
    }
    for (auto klist = 0u; klist < numNodeLists; ++klist) {
      neighborsijk.clear();
      if (domainDecompIndependent) {
        std::set_intersection(neighborsi[klist].begin(), neighborsi[klist].end(),
                              neighborsj[klist].begin(), neighborsj[klist].end(),
                              std::back_inserter(neighborsijk),
                              [&](const int a, const int b) { return mKeys(klist, a) < mKeys(klist, b); });
      } else {
        std::set_intersection(neighborsi[klist].begin(), neighborsi[klist].end(),
                              neighborsj[klist].begin(), neighborsj[klist].end(),
                              std::back_inserter(neighborsijk));
      }
      if (usePosition) {
        std::copy_if(neighborsijk.begin(), neighborsijk.end(), std::back_inserter(result[klist]),
                     [&](int k) { return (closestPointOnSegment(position(klist, k), posi, posj, b)); });
      } else {
        result[klist] = neighborsijk;
      }
    }
  } else if (i < (int)firstGhostNodei) {
    result = this->connectivityForNode(nodeListi, i);
  } else {
    result = this->connectivityForNode(nodeListj, j);
  }
  result[nodeListi].push_back(i);
  result[nodeListj].push_back(j);

  // That's it.
  TIME_END("ConnectivityMap_computeIntersectionConnectivity");
  return result;
}

//------------------------------------------------------------------------------
// Compute the union of neighbors for a pair of nodes.  Note this method 
// returns by value since this information is not stored by ConnectivityMap.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<vector<int> >
ConnectivityMap<Dimension>::
connectivityUnionForNodes(const int nodeListi, const int i,
                          const int nodeListj, const int j) const {

  // Pre-conditions.
  const unsigned numNodeLists = mNodeLists.size();
  const auto domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
  const auto ghostConnectivity = (mBuildGhostConnectivity or
                                  mBuildOverlapConnectivity or
                                  domainDecompIndependent);
  const unsigned firstGhostNodei = mNodeLists[nodeListi]->firstGhostNode();
  const unsigned firstGhostNodej = mNodeLists[nodeListj]->firstGhostNode();
  CONTRACT_VAR(ghostConnectivity);
  CONTRACT_VAR(firstGhostNodei);
  CONTRACT_VAR(firstGhostNodej);
  REQUIRE(nodeListi < (int)numNodeLists and
          nodeListj < (int)numNodeLists);
  REQUIRE(ghostConnectivity or i < (int)firstGhostNodei or j < (int)firstGhostNodej);

  // Do the deed.
  vector<vector<int> > result(numNodeLists);
  vector<vector<int> > neighborsi = this->connectivityForNode(nodeListi, i);
  vector<vector<int> > neighborsj = this->connectivityForNode(nodeListj, j);
  CHECK(neighborsi.size() == numNodeLists);
  CHECK(neighborsj.size() == numNodeLists);
  for (unsigned k = 0; k != numNodeLists; ++k) {
    sort(neighborsi[k].begin(), neighborsi[k].end());
    sort(neighborsj[k].begin(), neighborsj[k].end());
    set_union(neighborsi[k].begin(), neighborsi[k].end(),
              neighborsj[k].begin(), neighborsj[k].end(),
              back_inserter(result[k]));
  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Return the connectivity in terms of global node IDs.
//------------------------------------------------------------------------------
template<typename Dimension>
map<int, vector<int> > 
ConnectivityMap<Dimension>::
globalConnectivity(vector<Boundary<Dimension>*>& boundaries) const {

  // Get the set of global node IDs.
  auto globalIDs = globalNodeIDs<Dimension, typename vector<const NodeList<Dimension>*>::const_iterator>
    (mNodeLists.begin(), mNodeLists.end());

  // Make sure all ghost nodes have the appropriate global IDs.
  for (typename vector<Boundary<Dimension>*>::iterator boundItr = boundaries.begin();
       boundItr != boundaries.end();
       ++boundItr) (*boundItr)->applyFieldListGhostBoundary(globalIDs);
  for (typename vector<Boundary<Dimension>*>::iterator boundItr = boundaries.begin();
       boundItr != boundaries.end();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Now convert our connectivity to global IDs.
  map<int, vector<int> > result;
  const size_t numNodeLists = mNodeLists.size();
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {

    const NodeList<Dimension>* nodeListPtr = mNodeLists[nodeListi];
    for (auto i = 0u; i != nodeListPtr->numInternalNodes(); ++i) {
      const int gid = globalIDs(nodeListi, i);
      result[gid] = vector<int>();

      const vector< vector<int> >& fullConnectivity = connectivityForNode(nodeListPtr, i);
      CHECK(fullConnectivity.size() == numNodeLists);
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& connectivity = fullConnectivity[nodeListj];

        for (typename vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) result[gid].push_back(globalIDs(nodeListj, *jItr));

      }
      ENSURE(result[gid].size() == numNeighborsForNode(nodeListPtr, i));
    }
  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Compute the index for the given NodeList in our known set.
//------------------------------------------------------------------------------
template<typename Dimension>
unsigned
ConnectivityMap<Dimension>::
nodeListIndex(const NodeList<Dimension>* nodeListPtr) const {
  return distance(mNodeLists.begin(), 
                  find(mNodeLists.begin(), mNodeLists.end(), nodeListPtr));
}

//------------------------------------------------------------------------------
// Valid test.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ConnectivityMap<Dimension>::
valid() const {
  TIME_BEGIN("ConnectivityMap_valid");

  const auto domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
  const auto ghostConnectivity = (mBuildGhostConnectivity or
                                  mBuildOverlapConnectivity or
                                  domainDecompIndependent);
  const auto numNodeLists = mNodeLists.size();

  // We have to have the pairs
  if (not mNodePairListPtr) {
    cerr << "ConnectivityMap::valid: failed pairs existence check" << endl;
    return false;
  }

  // Make sure that the NodeLists are listed in the correct sequence, and are
  // FluidNodeLists.
  {
    const auto& registrar = NodeListRegistrar<Dimension>::instance();
    const auto  names = registrar.registeredNames();
    int lastPosition = -1;
    for (auto* nptr: mNodeLists) {
      const int newPosition = distance(names.begin(),
                                       find(names.begin(), names.end(), nptr->name()));
      if (newPosition <= lastPosition) {
        cerr << "ConnectivityMap::valid: Failed ordering of NodeLists at " << nptr->name() << endl;
        return false;
      }
      lastPosition = newPosition;
    }
  }

  // The following checks only apply to the per-point connectivity if available
  if (not mConnectivity.empty()) {

    // Check the offsets.
    if (mOffsets.size() != numNodeLists) {
      cerr << "ConnectivityMap::valid: Failed mOffsets.size() == numNodeLists" << endl;
      return false;
    }
    {
      const auto numNodes = (ghostConnectivity ? 
                             mNodeLists.back()->numNodes() : 
                             mNodeLists.back()->numInternalNodes());
      if (mConnectivity.size() != mOffsets.back() + numNodes) {
        cerr << "ConnectivityMap::valid: Failed offset bounding: " << mConnectivity.size() << " != " << mOffsets.back() << " + " << numNodes << endl;
        return false;
      }
    }

    // Iterate over each NodeList entered.
    for (auto [nodeListIDi, nodeListPtri]: enumerate(mNodeLists)) {

      // Are all internal nodes represented?
      const auto numNodes = (ghostConnectivity ?
                             nodeListPtri->numNodes() : 
                             nodeListPtri->numInternalNodes());
      //const int firstGhostNodei = nodeListPtri->firstGhostNode();
      if (((nodeListIDi < numNodeLists - 1u) and ((mOffsets[nodeListIDi + 1] - mOffsets[nodeListIDi]) != numNodes)) or
          ((nodeListIDi == numNodeLists - 1u) and ((mConnectivity.size() - (size_t)mOffsets[nodeListIDi]) != numNodes))) {
        cerr << "ConnectivityMap::valid: Failed test that all nodes set for NodeList "
             << mNodeLists[nodeListIDi]->name()
             << endl;
        return false;
      }

      // Iterate over the nodes for this NodeList.
      const int ioff = mOffsets[nodeListIDi];
      for (auto i = 0u; i < numNodes; ++i) {

        // The set of neighbors for this node.  This has to be sized as the number of
        // NodeLists.
        const vector< vector<int> >& allNeighborsForNode = mConnectivity[ioff + i];
        if (allNeighborsForNode.size() != numNodeLists) {
          cerr << "ConnectivityMap::valid: Failed allNeighborsForNode.size() == numNodeLists" << endl;
          return false;
        }

        // Iterate over the sets of NodeList neighbors for this node.
        for (auto nodeListIDj = 0u; nodeListIDj < numNodeLists; ++nodeListIDj) {
          const NodeList<Dimension>* nodeListPtrj = mNodeLists[nodeListIDj];
          //const int firstGhostNodej = nodeListPtrj->firstGhostNode();
          const vector<int>& neighbors = allNeighborsForNode[nodeListIDj];

          // We require that the node IDs be sorted, unique, and of course in a valid range.
          if (neighbors.size() > 0) {
            const auto minNeighbor = *min_element(neighbors.begin(), neighbors.end());
            const auto maxNeighbor = *max_element(neighbors.begin(), neighbors.end());

            if (minNeighbor < 0 or (size_t)maxNeighbor >= nodeListPtrj->numNodes()) {
              cerr << "ConnectivityMap::valid: Failed test that neighbors must be valid IDs: " << minNeighbor << " " << maxNeighbor << " " << nodeListPtrj->numNodes() << endl;
              return false;
            }

            // // When enforcing domain independence the ith node may be a ghost, but all of it's neighbors should
            // // be internal.
            // if (domainDecompIndependent and (i >= firstGhostNodei) and (maxNeighbor > firstGhostNodej)) {
            //   cerr << "ConnectivityMap::valid: Failed test that all neighbors of a ghost node should be internal." << endl;
            //   return false;
            // }

            for (auto k = 1u; k < neighbors.size(); ++k) {
              if (domainDecompIndependent) {
                // In the case of domain decomposition reproducibility, neighbors are sorted
                // by hashed IDs.
                if (mKeys(nodeListIDj, neighbors[k]) < mKeys(nodeListIDj, neighbors[k - 1])) {
                  cerr << "ConnectivityMap::valid: Failed test that neighbors must be sorted for node "
                       << i << endl;
                  for (vector<int>::const_iterator itr = neighbors.begin();
                       itr != neighbors.end();
                       ++itr) cerr << "(" << *itr << " " << mKeys(nodeListIDj, *itr) << ") ";
                  cerr << endl;
                  return false;
                }

              } else {
                // Otherwise they should be sorted by local ID.
                if (neighbors[k] <= neighbors[k - 1]) {
                  cerr << "ConnectivityMap::valid: Failed test that neighbors must be sorted" << endl;
                  for (vector<int>::const_iterator itr = neighbors.begin();
                       itr != neighbors.end();
                       ++itr) cerr << " " << *itr;
                  cerr << endl;
                  return false;
                }
              }
            }
          }

          // Check that the connectivity is symmetric.
          for (auto j: neighbors) {
            if (ghostConnectivity or ((size_t)j < nodeListPtrj->numInternalNodes())) {
              const vector< vector<int> >& otherNeighbors = connectivityForNode(nodeListPtrj, j);
              if (find(otherNeighbors[nodeListIDi].begin(),
                       otherNeighbors[nodeListIDi].end(),
                       i) == otherNeighbors[nodeListIDi].end()) {
                cerr << "ConnectivityMap::valid: Failed test that neighbors must be symmetric: " 
                     << i << " <> " << j 
                     << "  numneigbors(i)=" << neighbors.size() 
                     << "  numneigbors(j)=" << otherNeighbors[nodeListIDi].size() 
                     << endl;
                cerr << "   " << i << " : ";
                std::copy(neighbors.begin(), neighbors.end(), std::ostream_iterator<int>(std::cerr, " "));
                cerr << endl
                     << "   " << j << " : ";
                std::copy(otherNeighbors[nodeListIDi].begin(), otherNeighbors[nodeListIDi].end(), std::ostream_iterator<int>(std::cerr, " "));
                cerr << endl;
                return false;
              }
            }
          }

        }
      }
    }
  }

  // Check that if we are using domain decompostion independence then the keys 
  // have been calculated.
  if (domainDecompIndependent) {
    for (auto* nptr: mNodeLists) {
      if (not mKeys.haveNodeList(*nptr)) {
        cerr << "ConnectivityMap::valid: missing information from Keys for " <<  nptr->name() << endl;
        return false;
      }
    }
  }

  // Make sure all nodes are listed in the node index traversal stuff.
  if (mNodeTraversalIndices.size() != mNodeLists.size()) {
    cerr << "ConnectivityMap::valid: mNodeTraversalIndices wrong size!" << endl;
    return false;
  }
  for (auto nodeList = 0u; nodeList < numNodeLists; ++nodeList) {
    const auto numExpected = domainDecompIndependent ? mNodeLists[nodeList]->numNodes() : mNodeLists[nodeList]->numInternalNodes();
    bool ok = mNodeTraversalIndices[nodeList].size() == numExpected;
    for (auto i = 0u; i < numExpected; ++i) {
      ok = ok and (count(mNodeTraversalIndices[nodeList].begin(),
                         mNodeTraversalIndices[nodeList].end(),
                         i) == 1);
    }
    if (not ok) {
      cerr << "ConnectivityMap::valid: mNodeTraversalIndices elements messed up!" << endl;
      return false;
    }
  }

  // Everything must be OK.
  TIME_END("ConnectivityMap_valid");
  return true;
}

//------------------------------------------------------------------------------
// Internal method to build the connectivity for the requested set of NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConnectivityMap<Dimension>::
computeConnectivity() {
  TIME_BEGIN("ConnectivityMap_computeConnectivity");

  using Key = typename KeyTraits::Key;

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    for (auto* nptr: mNodeLists) {
      CONTRACT_VAR(nptr);
      REQUIRE(nptr->neighbor().valid());
    }
  }
  END_CONTRACT_SCOPE

  // Do we need to build the ghost connectivity as well?
  const bool domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
  const bool ghostConnectivity = (mBuildGhostConnectivity or
                                  mBuildOverlapConnectivity or
                                  domainDecompIndependent);

  // Build ourselves a temporary DataBase with the set of NodeLists.
  // Simultaneously find the maximum kernel extent.
  DataBase<Dimension> dataBase;
  auto kernelExtent = 0.0;
  for (auto* nptr: mNodeLists) {
    dataBase.appendNodeList(const_cast<NodeList<Dimension>&>(*nptr));
    kernelExtent = max(kernelExtent, (*nptr).neighbor().kernelExtent());
  }
  const auto kernelExtent2 = kernelExtent*kernelExtent;

  // Erase any prior information.
  const size_t numNodeLists = dataBase.numNodeLists();
  mOffsets.clear();
  mConnectivity.clear();
  mNodeTraversalIndices = vector<vector<int>>(numNodeLists);
  mIntersectionConnectivity.clear();

  // If we're trying to be domain decomposition independent, we need a key to sort
  // by that will give us a unique ordering regardless of position.  The Morton ordered
  // hash fills the bill.
  if (domainDecompIndependent) mKeys = mortonOrderIndices(dataBase);

  // Fill in the ordering for walking the nodes.
  if (domainDecompIndependent) {
    for (auto [k, nptr]: enumerate(mNodeLists)) {
      const auto n = nptr->numNodes();
      mNodeTraversalIndices[k].resize(n);
      vector<pair<int, Key>> keys;
      keys.reserve(n);
      for (auto i = 0u; i < n; ++i) keys.push_back(make_pair(i, mKeys(k, i)));
      sort(keys.begin(), keys.end(), ComparePairsBySecondElement<pair<int, Key> >());
      for (auto i = 0u; i < n; ++i) mNodeTraversalIndices[k][i] = keys[i].first;
      CHECK(mNodeTraversalIndices[k].size() == n);
    }
  } else {
    for (auto [k, nptr]: enumerate(mNodeLists)) {
      const auto n = nptr->numInternalNodes();
      mNodeTraversalIndices[k].resize(n);
      for (auto i = 0u; i < n; ++i) mNodeTraversalIndices[k][i] = i;
    }
  }

  // Create a list of flags to keep track of which nodes have been completed thus far.
  FieldList<Dimension, int> flagNodeDone = dataBase.newGlobalFieldList(int());
  flagNodeDone = 0;

  // Get the position and H fields.
  const auto position = dataBase.globalPosition();
  const auto H = dataBase.globalHfield();

  std::vector<NodePairIdxType> nodePairs;
  if (mNodePairListPtr) {
    nodePairs.reserve(mNodePairListPtr->size());
  }
  for (auto [iiNodeList, nptr]: enumerate(mNodeLists)) {
    const auto etaMax = nptr->neighbor().kernelExtent();

    // Iterate over the nodes in this NodeList, and look for any that are not done yet.
    const auto nii = (ghostConnectivity ?
                      nptr->numNodes() :
                      nptr->numInternalNodes());
    for (auto ii = 0u; ii < nii; ++ii) {
      if (flagNodeDone(iiNodeList, ii) == 0) {

        // Set the master nodes.
        vector<vector<int>> masterLists, coarseNeighbors;
        Neighbor<Dimension>::setMasterNeighborGroup(position(iiNodeList, ii),
                                                    H(iiNodeList, ii),
                                                    mNodeLists.begin(),
                                                    mNodeLists.end(),
                                                    etaMax,
                                                    masterLists,
                                                    coarseNeighbors,
                                                    ghostConnectivity);

        // Iterate over the full of NodeLists again to work on the master nodes.
        for (auto iNodeList = 0u; iNodeList != numNodeLists; ++iNodeList) {
          const auto nmaster = masterLists[iNodeList].size();
#pragma omp parallel 
          {
            std::vector<NodePairIdxType> nodePairs_private;
#pragma omp for schedule(dynamic)
            for (auto k = 0u; k < nmaster; ++k) {
              const auto i = masterLists[iNodeList][k];
              CHECK2(flagNodeDone(iNodeList, i) == 0, "(" << iNodeList << " " << i << ")");

              // Get the state for this node.
              const auto& ri = position(iNodeList, i);
              const auto& Hi = H(iNodeList, i);
              auto&       worki = mNodeLists[iNodeList]->work();
              const auto start = Timing::currentTime();

              // We keep track of the Morton indices.
              vector<vector<pair<int, Key>>> keys(numNodeLists);

              // Iterate over the neighbor NodeLists.
              for (auto jNodeList = 0u; jNodeList != numNodeLists; ++jNodeList) {
                const auto firstGhostNodej = mNodeLists[jNodeList]->firstGhostNode();

                // Iterate over the coarse neighbors in this NodeList.
                for (const auto j:  coarseNeighbors[jNodeList]) {
                  const auto& rj = position(jNodeList, j);
                  const auto& Hj = H(jNodeList, j);

                  // Compute the normalized distance between this pair.
                  const auto rij = ri - rj;
                  const auto eta2i = (Hi*rij).magnitude2();
                  const auto eta2j = (Hj*rij).magnitude2();

                  // If this pair is significant, add it to the list.
                  if (eta2i <= kernelExtent2 or eta2j <= kernelExtent2) {

                    // We don't include self-interactions.
                    if ((iNodeList != jNodeList) or (i != j)) {
                      if (calculatePairInteraction(iNodeList, i, jNodeList, j, firstGhostNodej)) nodePairs_private.push_back(NodePairIdxType(i, iNodeList, j, jNodeList));
                      if (domainDecompIndependent) keys[jNodeList].push_back(std::make_pair(j, mKeys(jNodeList, j)));
                    }
                  }
                }
              }
              CHECK(keys.size() == numNodeLists);
        
              // Flag this master node as done.
              flagNodeDone(iNodeList, i) = 1;
              worki(i) += Timing::difference(start, Timing::currentTime());
            }
            
            // Merge the NodePairList
#pragma omp critical
            nodePairs.insert(nodePairs.end(), nodePairs_private.begin(), nodePairs_private.end());
          } // end OMP parallel
        }
      }
    }
  }
  mNodePairListPtr = std::make_shared<NodePairList>(std::move(nodePairs));

  // Sort the NodePairList in order to enforce domain decomposition independence.
  if (domainDecompIndependent) {
    sortPairs(*mNodePairListPtr, mKeys);
  } else {
    std::sort(mNodePairListPtr->begin(), mNodePairListPtr->end());
  }

  // Overlap  and intersection connectivity are based on the per-point neighbors, so build that if necessary
  if (mBuildOverlapConnectivity or mBuildIntersectionConnectivity) buildPerPointConnectivity();

  // Do we need overlap connectivity?
  if (mBuildOverlapConnectivity) {
    TIME_BEGIN("ConnectivityMap_computeOverlapConnectivity");

    // To start out, *all* neighbors of a node (gather and scatter) are overlap neighbors.  Therefore we
    // first just copy the neighbor connectivity.
    mOverlapConnectivity = mConnectivity;

    for (auto iNodeList = 0u; iNodeList < numNodeLists; ++iNodeList) {
      const auto* nodeListPtr = mNodeLists[iNodeList];
      for (auto i = 0u; i < nodeListPtr->numNodes(); ++i) {
        const auto& neighborsi = mConnectivity[mOffsets[iNodeList] + i];
        CHECK(neighborsi.size() == numNodeLists);
        const auto& ri = position(iNodeList, i);
        const auto& Hi = H(iNodeList, i);

        // The points that have i in common are overlap neighbors with one another
        // for (auto jN1 = 0; jN1 < numNodeLists; ++jN1) {
        //   for (const auto j1 : neighborsi[jN1]) {
        //     for (auto jN2 = 0; jN2 < numNodeLists; ++jN2) {
        //       for (const auto j2 : neighborsi[jN2]) {
        //         if (!(jN1 == jN2  && j1 == j2)) {
        //           insertUnique(mOffsets, mOverlapConnectivity, mKeys, domainDecompIndependent,
        //                        jN1, j1, jN2, j2);
        //         }
        //       }
        //     }
        //   }
        // }
        
        // Find all the gather neighbors of i.
        for (auto jN1 = 0u; jN1 < numNodeLists; ++jN1) {
          for (const auto j1: neighborsi[jN1]) {
            const auto& rj1 = position(jN1, j1);
            const auto& Hj1 = H(jN1, j1);
            if ((Hi*(rj1 - ri)).magnitude2() <= kernelExtent2) {                           // Is j1 a gather neighbor of i?

              // Check if i and j1 have overlap directly.
              if ((Hj1*(rj1 - ri)).magnitude2() <= kernelExtent2) {
                insertUnique(mOffsets, mOverlapConnectivity, mKeys, domainDecompIndependent,
                             iNodeList, i, jN1, j1);
                insertUnique(mOffsets, mOverlapConnectivity, mKeys, domainDecompIndependent,
                             jN1, j1, iNodeList, i);
              }

              // Find the gather neighbors of j1, all of which share overlap with i.
              const auto& neighborsj1 = mConnectivity[mOffsets[jN1] + j1];
              for (auto jN2 = 0u; jN2 < numNodeLists; ++jN2) {
                for (const auto j2: neighborsj1[jN2]) {
                  const auto& rj2 = position(jN2, j2);
                  const auto& Hj2 = H(jN2, j2);
                  if ((Hj2*(rj2 - rj1)).magnitude2() <= kernelExtent2) {                   // Is j2 a scatter neighbor of j1?
                    insertUnique(mOffsets, mOverlapConnectivity, mKeys, domainDecompIndependent,
                                 iNodeList, i, jN2, j2);
                    insertUnique(mOffsets, mOverlapConnectivity, mKeys, domainDecompIndependent,
                                 jN2, j2, iNodeList, i);
                  }
                }
              }
            }
          }
        }
      }
    }
    TIME_END("ConnectivityMap_computeOverlapConnectivity");
  }

  // Are we building intersection connectivity?
  if (mBuildIntersectionConnectivity) {
    TIME_BEGIN("ConnectivityMap_precomputeIntersectionConnectivity");
    auto& pairs = *mNodePairListPtr;
    const auto npairs = pairs.size();
#pragma omp parallel
    {
      IntersectionConnectivityContainer intersection_private;
#pragma omp for
      for (auto k = 0u; k < npairs; ++k) {
        const auto& pair = pairs[k];
        intersection_private[pair] = this->connectivityIntersectionForNodes(pair.i_list, pair.i_node,
                                                                            pair.j_list, pair.j_node,
                                                                            position);
      }
#pragma omp critical
      {
        for (const auto& element: intersection_private) {
          mIntersectionConnectivity[element.first] = element.second;
        }
      } // omp critical
    }   // omp parallel
    TIME_END("ConnectivityMap_precomputeIntersectionConnectivity");
  }

  // {
  //   tpre = allReduce(unsigned(tpre), SPHERAL_OP_SUM) / Process::getTotalNumberOfProcesses() / CLOCKS_PER_SEC;
  //   tmaster = allReduce(unsigned(tmaster), SPHERAL_OP_SUM) / Process::getTotalNumberOfProcesses() / CLOCKS_PER_SEC;
  //   trefine = allReduce(unsigned(trefine), SPHERAL_OP_SUM) / Process::getTotalNumberOfProcesses() / CLOCKS_PER_SEC;
  //   twalk = allReduce(unsigned(twalk), SPHERAL_OP_SUM) / Process::getTotalNumberOfProcesses() / CLOCKS_PER_SEC;
  //   if (Process::getRank() == 0) {
  //     std::cerr << "ConnectivityMap timings (pre, master, refine, walk) = " << tpre << " " << tmaster << " " << trefine << " " << twalk << std::endl;
  //   }
  // }

  // Post conditions.
  BEGIN_CONTRACT_SCOPE
  // Make sure that the correct number of nodes have been completed.
  for (auto iNodeList = 0u; iNodeList < numNodeLists; ++iNodeList) {
    const auto n = (ghostConnectivity ? 
                    mNodeLists[iNodeList]->numNodes() :
                    mNodeLists[iNodeList]->numInternalNodes());
    for (auto i = 0u; i < n; ++i) {
      ENSURE2(flagNodeDone(iNodeList, i) == 1,
              "Missed connnectivity for (" << iNodeList << " " << i << ")");
    }
  }

  // Make sure we're ready to be used.
  ENSURE(valid());
  END_CONTRACT_SCOPE

  TIME_END("ConnectivityMap_computeConnectivity");
}

//------------------------------------------------------------------------------
// Build the per point connectivity when needed
// This method assumes the NodeListPair connectivity has already been
// constructed, and we just need to translate that information into the internal
// ConnectivityStorageType format (mConnectivity)
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConnectivityMap<Dimension>::
buildPerPointConnectivity() {
  TIME_BEGIN("ConnectivityMap_buildPerPointConnectivity");

  // Preconditions
  REQUIRE(mOffsets.empty());
  REQUIRE(mConnectivity.empty());

  // Count how many nodes we need to use per NodeList
  const bool domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
  const bool includeGhosts = domainDecompIndependent or mBuildGhostConnectivity or mBuildOverlapConnectivity;
  const auto numNodesInNodeList = [&](const size_t k) { return (includeGhosts ?
                                                                mNodeLists[k]->numNodes() :
                                                                mNodeLists[k]->numInternalNodes()); };
  const auto numNodeLists = mNodeLists.size();
  vector<size_t> numNodes(numNodeLists);
  for (auto k = 0u; k < numNodeLists; ++k) numNodes[k] = numNodesInNodeList(k);

  // Build the offsets
  mOffsets.resize(numNodeLists);
  mOffsets[0u] = 0u;
  for (auto k = 1u; k < numNodeLists; ++k) {
    mOffsets[k] = mOffsets[k - 1u] + numNodes[k - 1u];
  }

  // Build the empty per-point connectivity
  mConnectivity.resize(mOffsets.back() + numNodes.back(), vector<vector<int>>(numNodeLists));

  // Walk the pairs and build the per-point connectivity
  const auto& pairs = *mNodePairListPtr;
  for (const auto& pair: pairs) {
    const auto ki = pair.i_list;
    const auto kj = pair.j_list;
    const auto i = pair.i_node;
    const auto j = pair.j_node;
    CHECK(ki < numNodeLists);
    CHECK(kj < numNodeLists);
    CHECK(i < numNodesInNodeList(ki));
    CHECK(j < numNodesInNodeList(kj));

    CHECK(mOffsets[ki] + i < mConnectivity.size());
    auto& neighborsi = mConnectivity[mOffsets[ki] + i];
    CHECK(neighborsi.size() == numNodeLists);
    neighborsi[kj].push_back(j);

    if (j < mNodeLists[kj]->firstGhostNode() or includeGhosts) {
      CHECK(mOffsets[kj] + j < mConnectivity.size());
      auto& neighborsj = mConnectivity[mOffsets[kj] + j];
      CHECK(neighborsj.size() == numNodeLists);
      neighborsj[ki].push_back(i);
    }
  }

  TIME_END("ConnectivityMap_buildPerPointConnectivity");
}

}
