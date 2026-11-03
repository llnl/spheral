//---------------------------------Spheral++----------------------------------//
// ConnectivityMap
//
// Stores the full set of significant neighbors for a set of NodeLists.
//
// Created by J. Michael Owen, Sun Oct 30 15:36:33 PST 2005
//----------------------------------------------------------------------------//
#include "config.hh"
#include "ConnectivityMap.hh"
#include "NodeList/NodeList.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/mortonOrderIndices.hh"
#include "Utilities/PairComparisons.hh"
#include "Utilities/pointDistances.hh"
#include "Utilities/Timer.hh"
#include "Utilities/SpheralMessage.hh"

#include <ctime>

namespace Spheral {

namespace {
//------------------------------------------------------------------------------
// Present ConnectivityMap's owning host vectors through the common,
// non-owning NeighborViewContainer interface.
//------------------------------------------------------------------------------
template<typename Dimension, typename NeighborType>
inline
std::tuple<bool, size_t, size_t>
setMasterList(const NeighborType& neighbor,
              const typename Dimension::Vector& position,
              const typename Dimension::SymTensor& H,
              std::vector<std::pair<size_t, size_t>>& masterList,
              std::vector<std::pair<size_t, size_t>>& coarseNeighbors,
              const bool ghostConnectivity) {
  using NodeID = std::pair<size_t, size_t>;
#ifdef SPHERAL_UNIFIED_MEMORY
  NeighborViewContainer<NodeID> nativeMasterList(masterList.data(), masterList.size());
  NeighborViewContainer<NodeID> nativeCoarseNeighbors(coarseNeighbors.data(), coarseNeighbors.size());
#else
  NeighborViewContainer<NodeID> nativeMasterList = chai::makeManagedArray(masterList.data(), masterList.size(), chai::CPU, false);
  NeighborViewContainer<NodeID> nativeCoarseNeighbors = chai::makeManagedArray(coarseNeighbors.data(), coarseNeighbors.size(), chai::CPU, false);
#endif
  return neighbor.setMasterList(position, H,
                                nativeMasterList, nativeCoarseNeighbors,
                                ghostConnectivity);
}
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
             std::vector<std::vector<std::vector<size_t>>>& indices,
             const KeyContainer& keys,
             const bool useKeys,
             const size_t jN1, const size_t j1,
             const size_t jN2, const size_t j2) {
  if (jN1 != jN2 or j1 != j2) {
    auto& overlap = indices[offsets[jN1] + j1][jN2];
    std::vector<size_t>::iterator itr;
    if (useKeys) {
      itr = std::lower_bound(overlap.begin(), overlap.end(), j2,
                             [&](const auto a, const auto b) { return keys(jN2, a) < keys(jN2, b); });
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
// The input iterators must dereference to const NodeList*, and must be 
// const FluidNodeLists underneath.
//------------------------------------------------------------------------------
template<typename Dimension>
ConnectivityMap<Dimension>::
ConnectivityMap(const DataBase<Dimension>& dataBase,
                const bool buildGhostConnectivity,
                const bool buildOverlapConnectivity,
                const bool buildIntersectionConnectivity):
  mDataBase(dataBase),
  mBuildGhostConnectivity(buildGhostConnectivity or buildIntersectionConnectivity),
  mBuildOverlapConnectivity(buildOverlapConnectivity),
  mBuildIntersectionConnectivity(buildIntersectionConnectivity),
  mOffsets(),
  mConnectivity(),
  mNodeTraversalIndices(),
  mKeys(FieldStorageType::CopyFields),
  mCouplingPtr(std::make_shared<NodeCoupling>()) {

  // The private method does the grunt work of filling in the connectivity once we have
  // established the set of NodeLists.
  rebuild(mBuildGhostConnectivity, mBuildOverlapConnectivity, mBuildIntersectionConnectivity);
}

//------------------------------------------------------------------------------
// Rebuild for a given set of NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConnectivityMap<Dimension>::
rebuild(const bool buildGhostConnectivity,
        const bool buildOverlapConnectivity,
        const bool buildIntersectionConnectivity) {
  mBuildGhostConnectivity = buildGhostConnectivity or buildIntersectionConnectivity;
  mBuildOverlapConnectivity = buildOverlapConnectivity;
  mBuildIntersectionConnectivity = buildIntersectionConnectivity;
  mOffsets.clear();

  computeConnectivity();
  ENSURE(valid());
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
  const auto numNodeLists = mDataBase.numNodeLists();
  const auto& nodeLists = mDataBase.nodeListPtrs();
  if (domainDecompIndependent) {
    for (auto iNodeList = 0u; iNodeList < numNodeLists; ++iNodeList) {
      for (auto i = 0u; i < nodeLists[iNodeList]->numNodes(); ++i) {
        if (flags(iNodeList, i) == 0) mKeys(iNodeList, i) = KeyTraits::maxKey;
      }
    }
  }

  // Iterate over the traversal ordering and fix it
  for (auto iNodeList = 0u; iNodeList < numNodeLists; ++iNodeList) {
    const auto numNodes = ((domainDecompIndependent or mBuildGhostConnectivity or mBuildOverlapConnectivity) ? 
                           nodeLists[iNodeList]->numNodes() :
                           nodeLists[iNodeList]->numInternalNodes());

    std::vector<size_t> iNodesToKill;
    std::vector<std::pair<int, Key>> keys;
#pragma omp parallel
    {
      std::vector<size_t> iNodesToKill_thread;
      std::vector<std::pair<int, Key>> keys_thread;

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
  auto& pairs = *mNodePairListPtr;
  if (domainDecompIndependent) {
    sortPairs(pairs, mKeys);
  } else {
    std::sort(pairs.begin(), pairs.end());
  }

  // Rebuild the number of neighbors per point
  mNumNeighbors = 0u;
  for (const auto& p: pairs) {
    ++mNumNeighbors(p.i_list, p.i_node);
    ++mNumNeighbors(p.j_list, p.j_node);
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
        vector<vector<size_t>> newintersect(numNodeLists);
        for (auto klist = 0u; klist < numNodeLists; ++klist) {
          for (const auto k: oldintersect[klist]) {
            if (flags(klist, k) != 0) newintersect[klist].push_back(old2new(klist, k));
          }
        }
        intersection[pair] = newintersect;
      }
    }
    mIntersectionConnectivity = std::move(intersection);
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
std::vector<std::vector<size_t>>
ConnectivityMap<Dimension>::
connectivityIntersectionForNodes(const size_t nodeListi, const size_t i,
                                 const size_t nodeListj, const size_t j,
                                 const FieldList<Dimension, Vector>& position) const {

  // Pre-conditions.
  TIME_BEGIN("ConnectivityMap_computeIntersectionConnectivity");
  const auto& nodeLists = mDataBase.nodeListPtrs();
  const auto numNodeLists = nodeLists.size();
  const auto domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
  const auto ghostConnectivity = (mBuildGhostConnectivity or
                                  mBuildOverlapConnectivity or
                                  domainDecompIndependent);
  const auto firstGhostNodei = nodeLists[nodeListi]->firstGhostNode();
  const auto firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();
  const bool usePosition = (position.numFields() == numNodeLists);
  REQUIRE(nodeListi < numNodeLists and
          nodeListj < numNodeLists);
  REQUIRE(ghostConnectivity or i < firstGhostNodei or j < firstGhostNodej);
  REQUIRE(position.numFields() == numNodeLists or position.numFields() == 0);

  // Prepare the result.
  std::vector<std::vector<size_t>> result(numNodeLists);

  // If both nodes are internal, we simply intersect their neighbor lists.
  if (ghostConnectivity or (i < firstGhostNodei and j < firstGhostNodej)) {
    const auto& neighborsi = this->connectivityForNode(nodeListi, i);
    const auto& neighborsj = this->connectivityForNode(nodeListj, j);
    CHECK(neighborsi.size() == numNodeLists);
    CHECK(neighborsj.size() == numNodeLists);
    std::vector<size_t> neighborsijk;
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
                              [&](const auto a, const auto b) { return mKeys(klist, a) < mKeys(klist, b); });
      } else {
        std::set_intersection(neighborsi[klist].begin(), neighborsi[klist].end(),
                              neighborsj[klist].begin(), neighborsj[klist].end(),
                              std::back_inserter(neighborsijk));
      }
      if (usePosition) {
        std::copy_if(neighborsijk.begin(), neighborsijk.end(), std::back_inserter(result[klist]),
                     [&](auto k) { return (closestPointOnSegment(position(klist, k), posi, posj, b)); });
      } else {
        result[klist] = neighborsijk;
      }
    }
  } else if (i < firstGhostNodei) {
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
std::vector<std::vector<size_t>>
ConnectivityMap<Dimension>::
connectivityUnionForNodes(const size_t nodeListi, const size_t i,
                          const size_t nodeListj, const size_t j) const {

  // Pre-conditions.
  const auto& nodeLists = mDataBase.nodeListPtrs();
  const auto numNodeLists = nodeLists.size();
  const auto domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
  const auto ghostConnectivity = (mBuildGhostConnectivity or
                                  mBuildOverlapConnectivity or
                                  domainDecompIndependent);
  const auto firstGhostNodei = nodeLists[nodeListi]->firstGhostNode();
  const auto firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();
  CONTRACT_VAR(ghostConnectivity);
  CONTRACT_VAR(firstGhostNodei);
  CONTRACT_VAR(firstGhostNodej);
  REQUIRE(nodeListi < numNodeLists and
          nodeListj < numNodeLists);
  REQUIRE(ghostConnectivity or i < firstGhostNodei or j < firstGhostNodej);

  // Do the deed.
  std::vector<std::vector<size_t>> result(numNodeLists);
  std::vector<std::vector<size_t>> neighborsi = this->connectivityForNode(nodeListi, i);
  std::vector<std::vector<size_t>> neighborsj = this->connectivityForNode(nodeListj, j);
  CHECK(neighborsi.size() == numNodeLists);
  CHECK(neighborsj.size() == numNodeLists);
  for (auto k = 0u; k < numNodeLists; ++k) {
    std::sort(neighborsi[k].begin(), neighborsi[k].end());
    std::sort(neighborsj[k].begin(), neighborsj[k].end());
    std::set_union(neighborsi[k].begin(), neighborsi[k].end(),
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
map<size_t, vector<size_t>> 
ConnectivityMap<Dimension>::
globalConnectivity(std::vector<Boundary<Dimension>*>& boundaries) const {

  // Get the set of global node IDs.
  auto globalIDs = globalNodeIDs(mDataBase);

  // Make sure all ghost nodes have the appropriate global IDs.
  for (auto* bc: boundaries) bc->applyFieldListGhostBoundary(globalIDs);
  for (auto* bc: boundaries) bc->finalizeGhostBoundary();

  // Now convert our connectivity to global IDs.
  std::map<size_t, std::vector<size_t>> result;
  const auto& nodeLists = mDataBase.nodeListPtrs();
  const auto  numNodeLists = nodeLists.size();
  CONTRACT_VAR(numNodeLists);
  for (auto [nodeListi, nodeListPtr]: enumerate(nodeLists)) {
    for (auto i = 0u; i < nodeListPtr->numInternalNodes(); ++i) {
      const auto gid = globalIDs(nodeListi, i);
      result[gid].clear();
      const auto& fullConnectivity = connectivityForNode(nodeListPtr, i);
      CHECK(fullConnectivity.size() == numNodeLists);
      for (auto nodeListj = 0u; nodeListj < numNodeLists; ++nodeListj) {
        const auto& connectivity = fullConnectivity[nodeListj];
        for (const auto j: connectivity) result[gid].push_back(globalIDs(nodeListj, j));
      }
      ENSURE(result[gid].size() == numNeighborsForNode(nodeListPtr, i));
    }
  }

  // That's it.
  return result;
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
  const auto& nodeLists = mDataBase.nodeListPtrs();
  const auto numNodeLists = nodeLists.size();

  // We have to have the pairs
  if (not mNodePairListPtr) {
    SpheralWarning << "ConnectivityMap::valid: failed pairs existence check\n";
    return false;
  }

  // Make sure that the NodeLists are listed in the correct sequence, and are
  // FluidNodeLists.
  {
    const auto& registrar = NodeListRegistrar<Dimension>::instance();
    const auto  names = registrar.registeredNames();
    if (not std::ranges::is_sorted(names)) {
      SpheralWarning << "ConnectivityMap::valid: Failed ordering of NodeLists\n";
      return false;
    }
  }

  // The following checks only apply to the per-point connectivity if available
  if (not mConnectivity.empty()) {

    // Check the offsets.
    if (mOffsets.size() != numNodeLists) {
      SpheralWarning << "ConnectivityMap::valid: Failed mOffsets.size() == numNodeLists\n";
      return false;
    }
    {
      const auto numNodes = (ghostConnectivity ? 
                             nodeLists.back()->numNodes() : 
                             nodeLists.back()->numInternalNodes());
      if (mConnectivity.size() != mOffsets.back() + numNodes) {
        SpheralWarning << "ConnectivityMap::valid: Failed offset bounding: " << mConnectivity.size() << " != " << mOffsets.back() << " + " << numNodes << "\n";
        return false;
      }
    }

    // Iterate over each NodeList entered.
    for (auto [nodeListIDi, nodeListPtri]: enumerate(nodeLists)) {

      // Are all internal nodes represented?
      const auto numNodes = (ghostConnectivity ?
                             nodeListPtri->numNodes() : 
                             nodeListPtri->numInternalNodes());
      //const int firstGhostNodei = nodeListPtri->firstGhostNode();
      if (((nodeListIDi < numNodeLists - 1u) and ((mOffsets[nodeListIDi + 1] - mOffsets[nodeListIDi]) != numNodes)) or
          ((nodeListIDi == numNodeLists - 1u) and ((mConnectivity.size() - (size_t)mOffsets[nodeListIDi]) != numNodes))) {
        SpheralWarning << "ConnectivityMap::valid: Failed test that all nodes set for NodeList " << nodeListPtri->name() << "\n";
        return false;
      }

      // Iterate over the nodes for this NodeList.
      const auto ioff = mOffsets[nodeListIDi];
      for (auto i = 0u; i < numNodes; ++i) {

        // The set of neighbors for this node.  This has to be sized as the number of
        // NodeLists.
        const auto& allNeighborsForNode = mConnectivity[ioff + i];
        if (allNeighborsForNode.size() != numNodeLists) {
          SpheralWarning << "ConnectivityMap::valid: Failed allNeighborsForNode.size() == numNodeLists\n";
          return false;
        }

        // Iterate over the sets of NodeList neighbors for this node.
        for (auto [nodeListIDj, nodeListPtrj]: enumerate(nodeLists)) {
          //const int firstGhostNodej = nodeListPtrj->firstGhostNode();
          const auto& neighbors = allNeighborsForNode[nodeListIDj];

          // We require that the node IDs be sorted, unique, and of course in a valid range.
          if (neighbors.size() > 0) {
            const auto maxNeighbor = *std::max_element(neighbors.begin(), neighbors.end());

            if (maxNeighbor >= nodeListPtrj->numNodes()) {
              SpheralWarning << "ConnectivityMap::valid: Failed test that neighbors must be valid IDs: " << maxNeighbor << " " << nodeListPtrj->numNodes() << "\n";
              return false;
            }

            // // When enforcing domain independence the ith node may be a ghost, but all of it's neighbors should
            // // be internal.
            // if (domainDecompIndependent and (i >= firstGhostNodei) and (maxNeighbor > firstGhostNodej)) {
            //   SpheralWarning << "ConnectivityMap::valid: Failed test that all neighbors of a ghost node should be internal." << "\n";
            //   return false;
            // }

            for (auto k = 1u; k < neighbors.size(); ++k) {
              if (domainDecompIndependent) {
                // In the case of domain decomposition reproducibility, neighbors are sorted
                // by hashed IDs.
                if (mKeys(nodeListIDj, neighbors[k]) < mKeys(nodeListIDj, neighbors[k - 1])) {
                  SpheralWarning << "ConnectivityMap::valid: Failed test that neighbors must be sorted for node " << i << "\n";
                  for (const auto n: neighbors) SpheralWarning << "(" << n << " " << mKeys(nodeListIDj, n) << ") ";
                  SpheralWarning << "\n";
                  return false;
                }

              } else {
                // Otherwise they should be sorted by local ID.
                if (neighbors[k] <= neighbors[k - 1]) {
                  SpheralWarning << "ConnectivityMap::valid: Failed test that neighbors must be sorted" << "\n";
                  for (const auto n: neighbors) SpheralWarning << " " << n;
                  SpheralWarning << "\n";
                  return false;
                }
              }
            }
          }

          // Check that the connectivity is symmetric.
          for (auto j: neighbors) {
            if (ghostConnectivity or (j < nodeListPtrj->numInternalNodes())) {
              const auto& otherNeighbors = connectivityForNode(nodeListPtrj, j);
              if (std::find(otherNeighbors[nodeListIDi].begin(),
                            otherNeighbors[nodeListIDi].end(),
                            i) == otherNeighbors[nodeListIDi].end()) {
                SpheralWarning << "ConnectivityMap::valid: Failed test that neighbors must be symmetric: " 
                               << i << " <> " << j 
                               << "  numneigbors(i)=" << neighbors.size() 
                               << "  numneigbors(j)=" << otherNeighbors[nodeListIDi].size() 
                               << "\n";
                SpheralWarning << "   " << i << " : ";
                for (const auto n: neighbors) SpheralWarning << " " << n;
                SpheralWarning << "\n"
                               << "   " << j << " : ";
                for (const auto n: otherNeighbors[nodeListIDi]) SpheralWarning << " " << n;
                SpheralWarning << "\n";
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
    for (auto* nptr: nodeLists) {
      if (not mKeys.haveNodeList(*nptr)) {
        SpheralWarning << "ConnectivityMap::valid: missing information from Keys for " <<  nptr->name() << "\n";
        return false;
      }
    }
  }

  // Make sure all nodes are listed in the node index traversal stuff.
  if (mNodeTraversalIndices.size() != numNodeLists) {
    SpheralWarning << "ConnectivityMap::valid: mNodeTraversalIndices wrong size!" << "\n";
    return false;
  }
  for (auto nodeList = 0u; nodeList < numNodeLists; ++nodeList) {
    const auto numExpected = domainDecompIndependent ? nodeLists[nodeList]->numNodes() : nodeLists[nodeList]->numInternalNodes();
    bool ok = mNodeTraversalIndices[nodeList].size() == numExpected;
    for (auto i = 0u; i < numExpected; ++i) {
      ok = ok and (std::count(mNodeTraversalIndices[nodeList].begin(),
                              mNodeTraversalIndices[nodeList].end(),
                              i) == 1);
    }
    if (not ok) {
      SpheralWarning << "ConnectivityMap::valid: mNodeTraversalIndices elements messed up!" << "\n";
      return false;
    }
  }

  // // Check that the number of neighbors is set correctly
  // FieldList<Dimension, size_t> numNeighborsCheck(FieldStorageType::CopyFields);
  // for (const auto* nptr: mNodeLists) numNeighborsCheck.appendNewField("Number of neighbors", *nptr, 0u);
  // for (const auto& p: *mNodePairListPtr) {
  //   ++numNeighborsCheck(p.i_list, p.i_node);
  //   ++numNeighborsCheck(p.j_list, p.j_node);
  // }
  // if (mNumNeighbors != numNeighborsCheck) {
  //   cerr << "ConnectivityMap::valid: numNeighbors not tallied correctly:" << endl;
  //   for (auto k = 0u; k < numNodeLists; ++k) {
  //     const auto n = mNodeLists[k]->numNodes();
  //     for (auto i = 0u; i < n; ++i) {
  //       if (mNumNeighbors(k,i) != numNeighborsCheck(k,i)) cerr << "  @(" << k << "," << i << ") : " << mNumNeighbors(k,i) << " != " << numNeighborsCheck(k,i) << endl;
  //     }
  //   }
  //   return false;
  // }

  // Everything must be OK.
  TIME_END("ConnectivityMap_valid");
  return true;
}

//------------------------------------------------------------------------------
// computeConnectivity dispatch method necessary since Neighbor is a std::variant
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConnectivityMap<Dimension>::
computeConnectivity() {
  if (mDataBase.numNodeLists() == 0u) {
    mOffsets.clear();
    mConnectivity.clear();
    mOverlapConnectivity.clear();
    mNodeTraversalIndices.clear();
    mIntersectionConnectivity.clear();
    mNumNeighbors = mDataBase.newGlobalFieldList(size_t(0u), "Number of neighbors");
    mNodePairListPtr = std::make_shared<NodePairList>();
    return;
  }
  const auto& neighborVariant = mDataBase.neighbor();
  visitNeighbor(neighborVariant,
                [this](const auto& neighbor) {
                  this->computeConnectivityImpl(neighbor);
                });
}

//------------------------------------------------------------------------------
// Internal method to build the connectivity for the requested set of NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename NeighborType>
void
ConnectivityMap<Dimension>::
computeConnectivityImpl(const NeighborType& neighbor) {
  TIME_BEGIN("ConnectivityMap_computeConnectivity");

  using Key = typename KeyTraits::Key;

  // Do we need to build the ghost connectivity as well?
  const auto domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
  const auto ghostConnectivity = (mBuildGhostConnectivity or
                                  mBuildOverlapConnectivity or
                                  domainDecompIndependent);
  const auto nodeLists = mDataBase.nodeListPtrs();
  const auto numNodeLists = nodeLists.size();
  const auto etaMax = neighbor.kernelExtent();
  const auto etaMax2 = etaMax*etaMax;

  // Erase any prior information.
  mOffsets.clear();
  mConnectivity.clear();
  mOverlapConnectivity.clear();
  mNodeTraversalIndices = std::vector<std::vector<size_t>>(numNodeLists);
  mIntersectionConnectivity.clear();
  mNumNeighbors = mDataBase.newGlobalFieldList(size_t(0u), "Number of neighbors");

  // If we're trying to be domain decomposition independent, we need a key to sort
  // by that will give us a unique ordering regardless of position.  The Morton ordered
  // hash fills the bill.
  if (domainDecompIndependent) mKeys = mortonOrderIndices(mDataBase);

  // Fill in the ordering for walking the nodes.
  if (domainDecompIndependent) {
    for (auto [k, nptr]: enumerate(nodeLists)) {
      const auto n = nptr->numNodes();
      mNodeTraversalIndices[k].resize(n);
      std::vector<std::pair<size_t, Key>> keys;
      keys.reserve(n);
      for (size_t i = 0u; i < n; ++i) keys.push_back(make_pair(i, mKeys(k, i)));
      std::sort(keys.begin(), keys.end(), ComparePairsBySecondElement<std::pair<int, Key> >());
      for (size_t i = 0u; i < n; ++i) mNodeTraversalIndices[k][i] = keys[i].first;
      CHECK(mNodeTraversalIndices[k].size() == n);
    }
  } else {
    for (auto [k, nptr]: enumerate(nodeLists)) {
      const auto n = nptr->numInternalNodes();
      mNodeTraversalIndices[k].resize(n);
      for (size_t i = 0u; i < n; ++i) mNodeTraversalIndices[k][i] = i;
    }
  }

  // Create a list of flags to keep track of which nodes have been completed thus far.
  FieldList<Dimension, int> flagNodeDone = mDataBase.newGlobalFieldList(0);

  // Get the position and H fields.
  const auto position = mDataBase.globalPosition();
  const auto H = mDataBase.globalHfield();

  // Prepare memory for our results
  std::vector<NodePairIdxType> nodePairs;
  if (mNodePairListPtr) nodePairs.reserve(mNodePairListPtr->size());
  std::vector<NodeID> masterList(Dimension::pownu(8u)), coarseNeighbors(Dimension::pownu(20u));
  
  // Walk the NodeLists
  for (auto [iiNodeList, nptr]: enumerate(nodeLists)) {

    // Iterate over the nodes in this NodeList, and look for any that are not done yet.
    const auto nii = (ghostConnectivity ?
                      nptr->numNodes() :
                      nptr->numInternalNodes());
    for (auto ii = 0u; ii < nii; ++ii) {
      if (flagNodeDone(iiNodeList, ii) == 0) {

        // Set the master nodes.
        bool ok = false;
        size_t numMaster, numCoarse;
        while (not ok) {
          std::tie(ok, numMaster, numCoarse) = setMasterList<Dimension>(neighbor,
                                                                        position(iiNodeList, ii),
                                                                        H(iiNodeList, ii),
                                                                        masterList,
                                                                        coarseNeighbors,
                                                                        ghostConnectivity);
          if (masterList.size() < size_t(1.1*numMaster)) masterList.resize(size_t(1.1*numMaster));
          if (coarseNeighbors.size() < size_t(1.1*numCoarse)) coarseNeighbors.resize(size_t(1.1*numCoarse));
        }
        CHECK(ok);

        // Iterate over the master nodes
        for (auto imaster = 0u; imaster < numMaster; ++imaster) {
          const auto [iNodeList, i] = masterList[imaster];
          CHECK2(flagNodeDone(iNodeList, i) == 0, "(" << iNodeList << " " << i << ")");

          // Get the state for this node.
          const auto& ri = position(iNodeList, i);
          const auto& Hi = H(iNodeList, i);
          auto&       worki = nodeLists[iNodeList]->work();
          const auto start = Timing::currentTime();

          // Check the coarse neighbors
          for (auto icoarse = 0u; icoarse < numCoarse; ++icoarse) {
            const auto [jNodeList, j] = coarseNeighbors[icoarse];
            const auto firstGhostNodej = nodeLists[jNodeList]->firstGhostNode();
            const auto& rj = position(jNodeList, j);
            const auto& Hj = H(jNodeList, j);

            // Compute the normalized distance between this pair.
            const auto rij = ri - rj;
            const auto eta2i = (Hi*rij).magnitude2();
            const auto eta2j = (Hj*rij).magnitude2();

            // If this pair is significant, add it to the list.
            if (eta2i <= etaMax2 or eta2j <= etaMax2) {

              // We don't include self-interactions.
              if ((iNodeList != jNodeList or i != j) and
                  calculatePairInteraction(iNodeList, i, jNodeList, j, firstGhostNodej)) nodePairs.push_back(NodePairIdxType(i, iNodeList, j, jNodeList));
            }
          }
        
          // Flag this master node as done.
          flagNodeDone(iNodeList, i) = 1;
          worki(i) += Timing::difference(start, Timing::currentTime());
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

  // Build the number of neighbors per point
  for (const auto& p: *mNodePairListPtr) {
    ++mNumNeighbors(p.i_list, p.i_node);
    ++mNumNeighbors(p.j_list, p.j_node);
  }

  // Overlap  and intersection connectivity are based on the per-point neighbors, so build that if necessary
  if (mBuildOverlapConnectivity or mBuildIntersectionConnectivity) buildPerPointConnectivity();

  // Do we need overlap connectivity?
  if (mBuildOverlapConnectivity) {
    TIME_BEGIN("ConnectivityMap_computeOverlapConnectivity");

    // To start out, *all* neighbors of a node (gather and scatter) are overlap neighbors.  Therefore we
    // first just copy the neighbor connectivity.
    mOverlapConnectivity = mConnectivity;

    for (auto [iNodeList, nodeListPtr]: enumerate(nodeLists)) {
      const auto n = nodeListPtr->numNodes();
      for (auto i = 0u; i < n; ++i) {
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
            if ((Hi*(rj1 - ri)).magnitude2() <= etaMax2) {                           // Is j1 a gather neighbor of i?

              // Check if i and j1 have overlap directly.
              if ((Hj1*(rj1 - ri)).magnitude2() <= etaMax2) {
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
                  if ((Hj2*(rj2 - rj1)).magnitude2() <= etaMax2) {                   // Is j2 a scatter neighbor of j1?
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
  for (auto [iNodeList, nptr]: enumerate(nodeLists)) {
    const auto n = (ghostConnectivity ? 
                    nptr->numNodes() :
                    nptr->numInternalNodes());
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

  // cerr << "ConnectivityMap::buildPerPointConnectivity" << endl;

  // Double check we're starting fresh
  mOffsets.clear();
  mConnectivity.clear();

  // Count how many nodes we need to use per NodeList
  const bool domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
  const bool includeGhosts = domainDecompIndependent or mBuildGhostConnectivity or mBuildOverlapConnectivity;
  const auto& nodeLists = mDataBase.nodeListPtrs();
  const auto numNodesInNodeList = [&](const size_t k) { return (includeGhosts ?
                                                                nodeLists[k]->numNodes() :
                                                                nodeLists[k]->numInternalNodes()); };
  const auto numNodeLists = nodeLists.size();
  if (numNodeLists == 0u) {
    mOverlapConnectivity.clear();
    TIME_END("ConnectivityMap_buildPerPointConnectivity");
    return;
  }
  std::vector<size_t> numNodes(numNodeLists);
  for (auto k = 0u; k < numNodeLists; ++k) numNodes[k] = numNodesInNodeList(k);

  // Build the offsets
  mOffsets.resize(numNodeLists);
  mOffsets[0u] = 0u;
  for (auto k = 1u; k < numNodeLists; ++k) {
    mOffsets[k] = mOffsets[k - 1u] + numNodes[k - 1u];
  }

  // Build the empty per-point connectivity
  mConnectivity.resize(mOffsets.back() + numNodes.back(), std::vector<std::vector<size_t>>(numNodeLists));

  // Walk the pairs and build the per-point connectivity
  const auto& pairs = *mNodePairListPtr;
  for (const auto& pair: pairs) {
    const auto ki = pair.i_list;
    const auto kj = pair.j_list;
    const auto i = pair.i_node;
    const auto j = pair.j_node;
    CHECK(ki < numNodeLists);
    CHECK(kj < numNodeLists);

    CHECK(mOffsets[ki] + i < mConnectivity.size());
    auto& neighborsi = mConnectivity[mOffsets[ki] + i];
    CHECK(neighborsi.size() == numNodeLists);
    neighborsi[kj].push_back(j);

    if (j < nodeLists[kj]->firstGhostNode() or includeGhosts) {
      CHECK(mOffsets[kj] + j < mConnectivity.size());
      auto& neighborsj = mConnectivity[mOffsets[kj] + j];
      CHECK(neighborsj.size() == numNodeLists);
      neighborsj[ki].push_back(i);
    }
  }

  TIME_END("ConnectivityMap_buildPerPointConnectivity");
}

}
