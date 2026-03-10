//---------------------------------Spheral++----------------------------------//
// FlatConnectivity
//
// Creates global indices for each point and flattened connectivity indices
//----------------------------------------------------------------------------//
#include "FlatConnectivity.hh"

#include <algorithm>

#include "Boundary/ConstantBoundary.hh"
#include "Boundary/InflowOutflowBoundary.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "Geometry/CellFaceFlag.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Distributed/allReduce.hh"
#include "Utilities/DBC.hh"
#include "Utilities/Timer.hh"
#include "Utilities/globalNodeIDs.hh"

namespace Spheral {

namespace { // anonymous
template<typename Dimension>
int
numFluidNeighbors(const std::vector<std::vector<int>>& connectivity,
                  const DataBase<Dimension>& dataBase)
{
   auto numNeighbors = 0;
   auto nodeListj = 0;
   for (auto nodeListItrj = dataBase.fluidNodeListBegin();
        nodeListItrj != dataBase.fluidNodeListEnd();
        ++nodeListItrj, ++nodeListj) {
      numNeighbors += connectivity[nodeListj].size();
   }
   return numNeighbors;
}
} // end namespace anonymous

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
FlatConnectivity<Dimension>::
FlatConnectivity():
  mIndexingInitialized(false),
  mOverlapIndexingInitialized(false),
  mGlobalIndexingInitialized(false),
  mSurfaceIndexingInitialized(false),
  mBoundaryInformationInitialized(false),
  mNumLocalNodes(0),
  mNumInternalLocalNodes(0),
  mNumConnectivityNodes(0),
  mNumGlobalNodes(0),
  mNumBoundaryNodes(0),
  mNumGlobalBoundaryNodes(0),
  mFirstGlobalIndex(0),
  mLastGlobalIndex(0) {
}

//------------------------------------------------------------------------------
// Compute the indices
//------------------------------------------------------------------------------
template<typename Dimension>
void 
FlatConnectivity<Dimension>::
computeIndices(const DataBase<Dimension>& dataBase) {
  TIME_FUNCTION;
  VERIFY(fluidNodeListsFirst(dataBase));
  
  // Get information from DataBase
  const auto numNodeListsDB = dataBase.numFluidNodeLists();
  const auto numNodesDB = dataBase.numFluidNodes();
  const auto numInternalNodesDB = dataBase.numFluidInternalNodes();

  const auto& connectivity = dataBase.connectivityMap();
  const auto requireGhostConnectivity = connectivity.buildGhostConnectivity();

  // Store size information
  mNumLocalNodes = numNodesDB;
  mNumInternalLocalNodes = numInternalNodesDB;
  mNumConnectivityNodes = requireGhostConnectivity ? mNumLocalNodes : mNumInternalLocalNodes;
  
  // Get flattened local indices
  mNodeToLocalIndex.resize(numNodeListsDB);
  mLocalToNodeIndex.resize(mNumLocalNodes);
  {
    auto index = 0;
    auto nodeListi = 0;
    // Add the internal nodes to the flattened index
    for (auto nodeListItri = dataBase.fluidNodeListBegin();
         nodeListItri != dataBase.fluidNodeListEnd();
         ++nodeListItri, ++nodeListi) {
      const auto numNodesi = (*nodeListItri)->numNodes();
      const auto numInternalNodesi = (*nodeListItri)->numInternalNodes();
      mNodeToLocalIndex[nodeListi].resize(numNodesi);
      for (auto nodei = 0u; nodei < numInternalNodesi;
           ++nodei, ++index) {
        CHECK(index < mNumInternalLocalNodes);
        mNodeToLocalIndex[nodeListi][nodei] = index;
        mLocalToNodeIndex[index] = std::make_pair(nodeListi, nodei);
      }
    }
    // Add the boundary nodes to the flattened index
    nodeListi = 0;
    for (auto nodeListItri = dataBase.fluidNodeListBegin();
         nodeListItri != dataBase.fluidNodeListEnd();
         ++nodeListItri, ++nodeListi) {
      const auto numNodesi = (*nodeListItri)->numNodes();
      const auto numInternalNodesi = (*nodeListItri)->numInternalNodes();
      for (auto nodei = numInternalNodesi; nodei < numNodesi;
           ++nodei, ++index) {
        CHECK(mNumInternalLocalNodes <= index && index < mNumLocalNodes);
        mNodeToLocalIndex[nodeListi][nodei] = index;
        mLocalToNodeIndex[index] = std::make_pair(nodeListi, nodei);
      }
    }
    
    CHECK(index == mNumLocalNodes);
  }
  
  // Store the flattened connectivity
  mNumNeighbors.resize(mNumConnectivityNodes);
  mFlatToLocalIndex.resize(mNumConnectivityNodes);
  mLocalToFlatIndex.resize(mNumConnectivityNodes);
  {
    auto nodeListi = 0;
    for (auto nodeListItri = dataBase.fluidNodeListBegin();
         nodeListItri != dataBase.fluidNodeListEnd();
         ++nodeListItri, ++nodeListi) {
      const auto numNodesi = requireGhostConnectivity ? (*nodeListItri)->numNodes() : (*nodeListItri)->numInternalNodes();
      for (auto nodei = 0u; nodei < numNodesi; ++nodei) {
        // Get data from the connectivity map
        const auto connectivityi = connectivity.connectivityForNode(nodeListi, nodei);
        const auto locali = mNodeToLocalIndex[nodeListi][nodei];
        const auto numNeighborsi = numFluidNeighbors(connectivityi, dataBase);
        
        // Resize the arrays
        CHECK(locali < mNumConnectivityNodes);
        mNumNeighbors[locali] = numNeighborsi + 1;
        mFlatToLocalIndex[locali].resize(numNeighborsi + 1);
        mLocalToFlatIndex[locali].clear();
        mLocalToFlatIndex[locali].reserve(numNeighborsi + 1);
        
        // Add the point itself
        auto index = 0;
        mFlatToLocalIndex[locali][index] = locali;
        mLocalToFlatIndex[locali][locali] = index;
        ++index;

        // Add the other points
        auto nodeListj = 0;
        for (auto nodeListItrj = dataBase.fluidNodeListBegin();
             nodeListItrj != dataBase.fluidNodeListEnd();
             ++nodeListItrj, ++nodeListj) {
          for (auto nodej : connectivityi[nodeListj]) {
            const auto localj = mNodeToLocalIndex[nodeListj][nodej];
            CHECK(index < numNeighborsi + 1);
            mFlatToLocalIndex[locali][index] = localj;
            mLocalToFlatIndex[locali][localj] = index;
            ++index;
          }
        }
        
        CHECK(index == numNeighborsi + 1);
      }
    }
  }

  mIndexingInitialized = true;
  mGhostIndexingInitialized = requireGhostConnectivity;
  
  BEGIN_CONTRACT_SCOPE
  {
    // Make sure local and node indices are reversible
    for (auto locali = 0; locali < mNumLocalNodes; ++locali) {
      const auto pairi = localToNode(locali);
      const auto nodeListi = pairi.first;
      const auto nodei = pairi.second;
      CHECK(locali == nodeToLocal(nodeListi, nodei));
      CONTRACT_VAR(nodeListi);
      CONTRACT_VAR(nodei);
    }

    // Make sure that connectivity is reversible
    for (auto locali = 0; locali < mNumConnectivityNodes; ++locali) {
      // Make sure sizes match up (meaning we didn't get duplicate values in map)
      CHECK(size_t(locali) < mNumNeighbors.size() &&
            size_t(locali) < mFlatToLocalIndex.size() &&
            size_t(locali) < mLocalToFlatIndex.size());
      CHECK(mFlatToLocalIndex[locali].size() == size_t(mNumNeighbors[locali]) &&
            mLocalToFlatIndex[locali].size() == size_t(mNumNeighbors[locali]));
      
      const auto numNeighborsi = numNeighbors(locali);
      for (auto flatj = 0; flatj < numNeighborsi; ++flatj) {
        CHECK(flatj == localToFlat(locali, flatToLocal(locali, flatj)));
      }
    }
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Compute the flattened overlap indices
//------------------------------------------------------------------------------
template<typename Dimension>
void 
FlatConnectivity<Dimension>::
computeOverlapIndices(const DataBase<Dimension>& dataBase) {
  TIME_FUNCTION;
  VERIFY(mIndexingInitialized);
  
  // Get information from DataBase
  const auto numNodeListsDB = dataBase.numFluidNodeLists();
  const auto numNodesDB = dataBase.numFluidNodes();
  const auto numInternalNodesDB = dataBase.numFluidInternalNodes();
  CONTRACT_VAR(numNodeListsDB);
  CONTRACT_VAR(numNodesDB);
  CONTRACT_VAR(numInternalNodesDB);
  const auto& connectivity = dataBase.connectivityMap();
  const auto requireGhostConnectivity = connectivity.buildGhostConnectivity();
  VERIFY(connectivity.buildOverlapConnectivity());
  VERIFY(!requireGhostConnectivity || mGhostIndexingInitialized);
  
  // Make sure number of nodes has not changed since computing indices
  VERIFY(numNodesDB == size_t(mNumLocalNodes));
  VERIFY(numNodeListsDB == mNodeToLocalIndex.size());
  VERIFY(numInternalNodesDB == size_t(mNumInternalLocalNodes));
  
  // Store the flattened overlap connectivity
  mNumOverlapNeighbors.resize(mNumConnectivityNodes);
  mFlatOverlapToLocalIndex.resize(mNumConnectivityNodes);
  mLocalToFlatOverlapIndex.resize(mNumConnectivityNodes);
  {
    auto nodeListi = 0;
    for (auto nodeListItri = dataBase.fluidNodeListBegin();
         nodeListItri != dataBase.fluidNodeListEnd();
         ++nodeListItri, ++nodeListi) {
      const auto numNodesi = requireGhostConnectivity ? (*nodeListItri)->numNodes() : (*nodeListItri)->numInternalNodes();
      for (auto nodei = 0u; nodei < numNodesi; ++nodei) {
        // Get data from the connectivity map
        const auto connectivityi = connectivity.overlapConnectivityForNode(nodeListi, nodei);
        const auto locali = mNodeToLocalIndex[nodeListi][nodei];
        const auto numNeighborsi = numFluidNeighbors(connectivityi, dataBase);
        
        // Resize the arrays 
        mNumOverlapNeighbors[locali] = numNeighborsi + 1;
        mFlatOverlapToLocalIndex[locali].resize(numNeighborsi + 1);
        mLocalToFlatOverlapIndex[locali].clear();
        mLocalToFlatOverlapIndex[locali].reserve(numNeighborsi + 1);

        // Add the point itself
        auto index = 0;
        mFlatOverlapToLocalIndex[locali][index] = locali;
        mLocalToFlatOverlapIndex[locali][locali] = index;
        ++index;
        
        // Add the other points
        auto nodeListj = 0;
        for (auto nodeListItrj = dataBase.fluidNodeListBegin();
             nodeListItrj != dataBase.fluidNodeListEnd();
             ++nodeListItrj, ++nodeListj) {
          for (auto nodej : connectivityi[nodeListj]) {
            const auto localj = mNodeToLocalIndex[nodeListj][nodej];
            CHECK(index < numNeighborsi + 1);
            mFlatOverlapToLocalIndex[locali][index] = localj;
            mLocalToFlatOverlapIndex[locali][localj] = index;
            ++index;
          }
        }
      }
    }
  }
  
  mOverlapIndexingInitialized = true;

  BEGIN_CONTRACT_SCOPE
  {
    // Make sure that overlap connectivity is reversible
    // const auto numNodesToCheck = requireGhostConnectivity ? mNumLocalNodes : mNumInternalLocalNodes;
    for (auto locali = 0; locali < mNumConnectivityNodes; ++locali) {
      // Make sure sizes match up (meaning we didn't get duplicate values in map)
      CHECK(mFlatOverlapToLocalIndex[locali].size() == size_t(mNumOverlapNeighbors[locali]));
      CHECK(mLocalToFlatOverlapIndex[locali].size() == size_t(mNumOverlapNeighbors[locali]));
      
      const auto numNeighborsi = numOverlapNeighbors(locali);
      for (auto flatj = 0; flatj < numNeighborsi; ++flatj) {
        CHECK(flatj == localToFlatOverlap(locali, flatOverlapToLocal(locali, flatj)));
      }
    }
    
    // // Make sure all points that should be overlap neighbors are
    // const auto position = dataBase.fluidPosition();
    // const auto H = dataBase.fluidHfield();
    // for (auto locali = 0; locali < mNumLocalNodes; ++locali) {
    //   const auto numNeighborsi = numNeighbors(locali);
    //   for (auto flatj = 0; flatj < numNeighborsi; ++flatj) {
    //     const auto localj = flatToLocal(locali, flatj);
    //     for (auto flatk = 0; flatk < numNeighborsi; ++flatk) {
    //       const auto localk = flatToLocal(locali, flatk);
    //       // Is k a neighbor of j?
    //       if (mLocalToFlatOverlapIndex[localj].count(localk) < 1) {
    //         std::cout << "i\t" << locali << "\t";
    //         std::cout << "j\t" << localj << "\t";
    //         std::cout << "k\t" << localk << "\t";
    //         std::cout << "dist\t" << (position(0, localj) - position(0, localk)).magnitude();
    //         std::cout << std::endl;
    //       }
    //       // CHECK2(mLocalToFlatOverlapIndex[localj].count(localk) > 0); 
    //     }
    //   }
    // }
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Compute the global indices
//------------------------------------------------------------------------------
template<typename Dimension>
void 
FlatConnectivity<Dimension>::
computeGlobalIndices(const DataBase<Dimension>& dataBase,
                     const std::vector<Boundary<Dimension>*>& boundaries) {
  TIME_FUNCTION;
  VERIFY(mIndexingInitialized);
  
  // Get information from DataBase
  const auto numNodeListsDB = dataBase.numFluidNodeLists();
  const auto numNodesDB = dataBase.numFluidNodes();
  const auto numInternalNodesDB = dataBase.numFluidInternalNodes();
  const auto numGlobalNodesDB = dataBase.globalNumFluidInternalNodes();
  CONTRACT_VAR(numNodeListsDB);
  CONTRACT_VAR(numNodesDB);
  CONTRACT_VAR(numInternalNodesDB);
  CONTRACT_VAR(numGlobalNodesDB);

  // Make sure number of nodes has not changed since computing indices
  VERIFY(numNodesDB == size_t(mNumLocalNodes));
  VERIFY(numNodeListsDB == mNodeToLocalIndex.size());
  VERIFY(numInternalNodesDB == size_t(mNumInternalLocalNodes));
  
  // Get global indices manually
  int globalScan = distScan(mNumInternalLocalNodes, SPHERAL_OP_SUM);
  VERIFY(globalScan >= mNumInternalLocalNodes);
  mFirstGlobalIndex = globalScan - mNumInternalLocalNodes;
  mLastGlobalIndex = globalScan - 1;
  mNumGlobalNodes = allReduce(mNumInternalLocalNodes, SPHERAL_OP_SUM);
  VERIFY(mNumGlobalNodes >= mNumInternalLocalNodes);
  VERIFY(mNumGlobalNodes == int(numGlobalNodesDB));
  // std::cout << Process::getRank() << "\t" << mNumInternalLocalNodes << "\t" << mNumGlobalNodes << "\t" << mFirstGlobalIndex << "\t" << mLastGlobalIndex << std::endl;
  
  FieldList<Dimension, int> globalNodeIndices = dataBase.newFluidFieldList(0, "global node IDs");
  // globalNodeIndices = globalNodeIDs(dataBase);

  // Fill the global node IDs
  auto currentGlobalIndex = mFirstGlobalIndex;
  for (auto locali = 0; locali < mNumInternalLocalNodes; ++locali) {
    const auto pairi = mLocalToNodeIndex[locali];
    const auto nodeListi = pairi.first;
    const auto nodei = pairi.second;
    globalNodeIndices(nodeListi, nodei) = currentGlobalIndex;
    currentGlobalIndex += 1;
  }
  VERIFY(currentGlobalIndex == mLastGlobalIndex + 1);
  
  // Apply boundary condition to indices
  for (auto* boundary : boundaries) {
    boundary->applyFieldListGhostBoundary(globalNodeIndices);
  }

  // Finalize boundary condition
  for (auto* boundary : boundaries) {
    boundary->finalizeGhostBoundary();
  }
  
  // Fill in the global indices
  mLocalToGlobalIndex.resize(mNumLocalNodes);
  for (auto locali = 0; locali < mNumLocalNodes; ++locali) {
    const auto pairi = mLocalToNodeIndex[locali];
    const auto nodeListi = pairi.first;
    const auto nodei = pairi.second;
    mLocalToGlobalIndex[locali] = globalNodeIndices(nodeListi, nodei);
  }
  
  mGlobalIndexingInitialized = true;
  
  BEGIN_CONTRACT_SCOPE
  {
    // Make sure the global indices are contiguous on this processor
    if (mNumInternalLocalNodes > 0) {
      auto prevGlobalIndex = mLocalToGlobalIndex[0];
      CONTRACT_VAR(prevGlobalIndex);
      for (auto locali = 1; locali < mNumInternalLocalNodes; ++locali) {
        const auto currGlobalIndex = mLocalToGlobalIndex[locali];
        CHECK(prevGlobalIndex == currGlobalIndex - 1);
        prevGlobalIndex = currGlobalIndex;
      }
    }
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Compute the surface indices
//------------------------------------------------------------------------------
template<typename Dimension>
void 
FlatConnectivity<Dimension>::
computeSurfaceIndices(const DataBase<Dimension>& dataBase,
                      const State<Dimension>& state) {
  TIME_FUNCTION;
  VERIFY(mIndexingInitialized);
  VERIFY(mGhostIndexingInitialized); // Could consider editing to not require this
  
  // Get information from the DataBase and State
  const auto numNodeListsDB = dataBase.numFluidNodeLists();
  const auto numNodesDB = dataBase.numFluidNodes();
  const auto numInternalNodesDB = dataBase.numFluidInternalNodes();
  CONTRACT_VAR(numNodeListsDB);
  CONTRACT_VAR(numNodesDB);
  CONTRACT_VAR(numInternalNodesDB);
  const auto& connectivity = dataBase.connectivityMap();
  const auto cells = state.fields(HydroFieldNames::cells, FacetedVolume());
  const auto cellFaceFlags = state.fields(HydroFieldNames::cellFaceFlags, std::vector<CellFaceFlag>());
  VERIFY(cells.size() == numNodeListsDB
         &&cellFaceFlags.size() == numNodeListsDB);
#if REPLACEOVERLAP
  const auto H = dataBase.fluidHfield();
  const auto position = dataBase.fluidPosition();
  const auto extent = dataBase.maxKernelExtent();
#endif
  
  // Make sure number of nodes has not changed since computing indices
  VERIFY(numNodesDB == size_t(mNumLocalNodes));
  VERIFY(numNodeListsDB == mNodeToLocalIndex.size());
  VERIFY(numInternalNodesDB == size_t(mNumInternalLocalNodes));

  // Since we are doing a gather operation, we need to make sure to clear out old data first
  mSurfaceNormal.resize(mNumLocalNodes);
  mSurfaceFlatIndex.resize(mNumLocalNodes);
  mVoidSurfaces.resize(mNumLocalNodes);
  for (auto i = 0; i < mNumLocalNodes; ++i) {
    mSurfaceNormal[i].clear();
    mSurfaceFlatIndex[i].clear();
    mVoidSurfaces[i].clear();
  }
  
  // For each cell, for each surface for that cell, add the surface to the
  // points that see the cell
  {
    ArrayDim normalArray;
    auto nodeListi = 0;
    for (auto nodeListItri = dataBase.fluidNodeListBegin();
         nodeListItri != dataBase.fluidNodeListEnd();
         ++nodeListItri, ++nodeListi) {
      const auto numNodesi = (*nodeListItri)->numNodes();
      for (auto nodei = 0u; nodei < numNodesi; ++nodei) {
        const auto locali = mNodeToLocalIndex[nodeListi][nodei];
        const auto& flags = cellFaceFlags(nodeListi, nodei);
        const auto numFlags = flags.size();
        if (numFlags > 0) {
          // Get the connectivity and surface information for this cell
          const auto& connectivityi = connectivity.connectivityForNode(nodeListi, nodei);
          const auto& cell = cells(nodeListi, nodei);
          const auto& facets = cell.facets();
#if REPLACEOVERLAP
          const auto xi = position(nodeListi, nodei);
#endif
        
          // If the surface flag is a void, then add the normal to the indexing
          for (auto flag : flags) {
            if (flag.nodeListj == -1) { // opposite side is a void
              mVoidSurfaces[locali].push_back(flag.cellFace);
              const auto& facet = facets[flag.cellFace];
              const auto& normalArea = facet.normal();
              const auto normal = normalArea.unitVector();

              normalToArray(normal, normalArray);
              for (auto nodeListj = 0u; nodeListj < numNodeListsDB; ++nodeListj) {
                for (auto nodej : connectivityi[nodeListj]) {
#if REPLACEOVERLAP
                  const auto xj = position(nodeListj, nodej);
                  const auto Hj = H(nodeListj, nodej);

                  // See if xi is inside support of point j
                  if (2.0 * (Hj * (xi - xj)).magnitude() <= extent) {
                    // If the normal is added, the index is the current size of the surface normal
                    const auto localj = mNodeToLocalIndex[nodeListj][nodej];
                    const auto index = mSurfaceNormal[localj].size();
                    const auto added = mSurfaceFlatIndex[localj].emplace(normalArray, index);
                    if (added.second) {
                      // We avoid roundoff error by adding the normal directly instead of recreating it from the indices
                      mSurfaceNormal[localj].push_back(normal);
                    }
                  }
#else
                  // If the normal is added, the index is the current size of the surface normal
                  const auto localj = mNodeToLocalIndex[nodeListj][nodej];
                  const auto index = mSurfaceNormal[localj].size();
                  const auto added = mSurfaceFlatIndex[localj].emplace(normalArray, index);
                  if (added.second) {
                    // We avoid roundoff error by adding the normal directly instead of recreating it from the indices
                    mSurfaceNormal[localj].push_back(normal);
                  }
#endif
                }
              }
              
              // Add the self-contribution
              const auto index = mSurfaceNormal[locali].size();
              const auto added = mSurfaceFlatIndex[locali].emplace(normalArray, index);
              if (added.second) {
                mSurfaceNormal[locali].push_back(normal);
              }
            }
          }
        }
      }
    }
  }
  
  mSurfaceIndexingInitialized = true;

  BEGIN_CONTRACT_SCOPE
  {
    for (auto locali = 0; locali < mNumLocalNodes; ++locali) {
      // Make sure the sizes of the map and normals match up
      const auto numSurfacesi = numSurfaces(locali);
      CHECK(mSurfaceNormal[locali].size() == size_t(numSurfacesi));
      CHECK(mSurfaceFlatIndex[locali].size() == size_t(numSurfacesi));

      // Make sure the map returns the expected surface index
      for (auto flats = 0; flats < numSurfacesi; ++flats) {
        const auto& normals = normal(locali, flats);
        const auto calcs = surfaceIndex(locali, normals);
        CHECK(calcs == flats);
        CONTRACT_VAR(calcs);
      }

      // Make sure all key/value pairs line up to an existing normal
      for (auto& norms : mSurfaceFlatIndex[locali]) {
          const auto& array = norms.first;
          const auto flats = norms.second;
          CHECK(flats < numSurfacesi);
          const auto& normal2 = normal(locali, flats);
          ArrayDim array2;
          normalToArray(normal2, array2);
          CHECK(array == array2);
          CONTRACT_VAR(array);
      }
    }
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Compute the boundary information
//------------------------------------------------------------------------------
template<typename Dimension>
void 
FlatConnectivity<Dimension>::
computeBoundaryInformation(const DataBase<Dimension>& dataBase,
                           const std::vector<Boundary<Dimension>*>& boundaries) {
  VERIFY(mIndexingInitialized);

  // Get information from the dataBase
  const auto numNodeListsDB = dataBase.numFluidNodeLists();
  const auto numNodesDB = dataBase.numFluidNodes();
  const auto numInternalNodesDB = dataBase.numFluidInternalNodes();
  CONTRACT_VAR(numNodeListsDB);
  CONTRACT_VAR(numNodesDB);
  CONTRACT_VAR(numInternalNodesDB);

  // Make sure the sizes haven't changed since the indexing was initialized
  VERIFY(numNodesDB == size_t(mNumLocalNodes));
  VERIFY(numNodeListsDB == mNodeToLocalIndex.size());
  VERIFY(numInternalNodesDB == size_t(mNumInternalLocalNodes));
  
  // Initialize the arrays
  mConstantBoundaryNodes.clear();
  mIsConstantBoundaryNode.assign(mNumLocalNodes, false);
  
  // Get the constant boundary nodes
  for (auto* boundary : boundaries) {
    if (dynamic_cast<const ConstantBoundary<Dimension>*>(boundary) != nullptr
        || dynamic_cast<const InflowOutflowBoundary<Dimension>*>(boundary) != nullptr) {
      mNumBoundaryNodes += boundary->numGhostNodes();
      auto nodeListi = 0;
      for (auto nodeListItr = dataBase.fluidNodeListBegin();
           nodeListItr != dataBase.fluidNodeListEnd();
           ++nodeListItr, ++nodeListi) {
        const auto ghostNodes = boundary->ghostNodes(**nodeListItr);
        for (auto nodei : ghostNodes) {
          const auto locali = mNodeToLocalIndex[nodeListi][nodei];
          CHECK(locali < mNumLocalNodes);
          mConstantBoundaryNodes.push_back(locali);
          mIsConstantBoundaryNode[locali] = true;
        }
      }
    }
  }
  
  // Make sure the constant boundary nodes are unique
  std::sort(mConstantBoundaryNodes.begin(), mConstantBoundaryNodes.end());
  mConstantBoundaryNodes.erase(std::unique(mConstantBoundaryNodes.begin(), mConstantBoundaryNodes.end()), mConstantBoundaryNodes.end());
  mNumBoundaryNodes = mConstantBoundaryNodes.size();
  mNumGlobalBoundaryNodes = allReduce(mNumBoundaryNodes, SPHERAL_OP_SUM);

  // For each point, get the number of neighbors that are constant boundary nodes
  {
    mNumConstantBoundaryNeighbors.resize(mNumConnectivityNodes);
    for (auto locali = 0; locali < mNumConnectivityNodes; ++locali) {
      auto num = 0;
      CHECK(size_t(locali) < mFlatToLocalIndex.size());
      for (auto localj : mFlatToLocalIndex[locali]) {
        CHECK(size_t(localj) < mIsConstantBoundaryNode.size());
        if (mIsConstantBoundaryNode[localj]) {
          num += 1;
        }
      }
      mNumConstantBoundaryNeighbors[locali] = num;
    }
  }

  // Do the same for overlap connectivity, if it has been computed
  if (mOverlapIndexingInitialized) {
    mNumConstantBoundaryOverlapNeighbors.resize(mNumConnectivityNodes);
    for (auto locali = 0; locali < mNumConnectivityNodes; ++locali) {
      auto num = 0;
      CHECK(size_t(locali) < mFlatOverlapToLocalIndex.size());
      for (auto localj : mFlatOverlapToLocalIndex[locali]) {
        CHECK(size_t(localj) < mIsConstantBoundaryNode.size());
        if (mIsConstantBoundaryNode[localj]) {
          num += 1;
        }
      }
      mNumConstantBoundaryOverlapNeighbors[locali] = num;
    }
  }

  mBoundaryInformationInitialized = true;
}

//------------------------------------------------------------------------------
// Get the local indices for the neighbors of point i
//------------------------------------------------------------------------------
template<typename Dimension>
void 
FlatConnectivity<Dimension>::
neighborIndices(const int locali,
                std::vector<int>& localNeighbors) const {
  CHECK(mIndexingInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(size_t(locali) < mFlatToLocalIndex.size());
  localNeighbors = mFlatToLocalIndex[locali];
}

template<typename Dimension>
void 
FlatConnectivity<Dimension>::
overlapNeighborIndices(const int locali,
                       std::vector<int>& localNeighbors) const {
  CHECK(mOverlapIndexingInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(size_t(locali) < mFlatOverlapToLocalIndex.size());
  localNeighbors = mFlatOverlapToLocalIndex[locali];
}

template<typename Dimension>
void 
FlatConnectivity<Dimension>::
constNeighborIndices(const int locali,
                std::vector<int>& localNeighbors) const {
  CHECK(mIndexingInitialized);
  CHECK(mBoundaryInformationInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(size_t(locali) < mFlatToLocalIndex.size());
  
  const auto numConstNeighbors = mNumConstantBoundaryNeighbors[locali];
  localNeighbors.resize(numConstNeighbors);
  auto index = 0;
  for (auto localj : mFlatToLocalIndex[locali]) {
    CHECK(size_t(localj) < mIsConstantBoundaryNode.size());
    if (mIsConstantBoundaryNode[localj]) {
      CHECK(index < numConstNeighbors);
      localNeighbors[index] = localj;
      ++index;
    }
  }
  CHECK(index == numConstNeighbors);
}

template<typename Dimension>
void 
FlatConnectivity<Dimension>::
overlapConstNeighborIndices(const int locali,
                       std::vector<int>& localNeighbors) const {
  CHECK(mOverlapIndexingInitialized);
  CHECK(mBoundaryInformationInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(size_t(locali) < mFlatOverlapToLocalIndex.size());

  const auto numConstNeighbors = mNumConstantBoundaryOverlapNeighbors[locali];
  localNeighbors.resize(numConstNeighbors);
  auto index = 0;
  for (auto localj : mFlatOverlapToLocalIndex[locali]) {
    CHECK(size_t(localj) < mIsConstantBoundaryNode.size());
    if (mIsConstantBoundaryNode[localj]) {
      CHECK(index < numConstNeighbors);
      localNeighbors[index] = localj;
      ++index;
    }
  }
  CHECK(index == numConstNeighbors);
}

template<typename Dimension>
void 
FlatConnectivity<Dimension>::
nonConstNeighborIndices(const int locali,
                        std::vector<int>& localNeighbors) const {
  CHECK(mIndexingInitialized);
  CHECK(mBoundaryInformationInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(size_t(locali) < mNumNeighbors.size() &&
        size_t(locali) < mFlatToLocalIndex.size() &&
        size_t(locali) < mNumConstantBoundaryNeighbors.size());

  const auto numNonConstNeighbors = mNumNeighbors[locali] - mNumConstantBoundaryNeighbors[locali];
  localNeighbors.resize(numNonConstNeighbors);
  auto index = 0;
  for (auto localj : mFlatToLocalIndex[locali]) {
    CHECK(size_t(localj) < mIsConstantBoundaryNode.size());
    if (!mIsConstantBoundaryNode[localj]) {
      CHECK(index < numNonConstNeighbors);
      localNeighbors[index] = localj;
      ++index;
    }
  }
  CHECK(index == numNonConstNeighbors);
}

template<typename Dimension>
void 
FlatConnectivity<Dimension>::
overlapNonConstNeighborIndices(const int locali,
                               std::vector<int>& localNeighbors) const {
  CHECK(mOverlapIndexingInitialized);
  CHECK(mBoundaryInformationInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(size_t(locali) < mFlatOverlapToLocalIndex.size());

  const auto numNonConstNeighbors = mNumOverlapNeighbors[locali] - mNumConstantBoundaryOverlapNeighbors[locali];
  localNeighbors.resize(numNonConstNeighbors);
  auto index = 0;
  for (auto localj : mFlatOverlapToLocalIndex[locali]) {
    CHECK(size_t(localj) < mIsConstantBoundaryNode.size());
    if (!mIsConstantBoundaryNode[localj]) {
      CHECK(index < numNonConstNeighbors);
      localNeighbors[index] = localj;
      ++index;
    }
  }
  CHECK(index == numNonConstNeighbors);
}

//------------------------------------------------------------------------------
// Get the global indices, without any constant points, for the point i
//------------------------------------------------------------------------------
template<typename Dimension>
void 
FlatConnectivity<Dimension>::
globalNeighborIndices(const int locali,
                      std::vector<int>& globalNeighbors) const {
  CHECK(mIndexingInitialized);
  CHECK(mBoundaryInformationInitialized);
  CHECK(mGlobalIndexingInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(size_t(locali) < mNumNeighbors.size() &&
        size_t(locali) < mFlatToLocalIndex.size());

  const auto numNonConstNeighbors = mNumNeighbors[locali] - mNumConstantBoundaryNeighbors[locali];
  globalNeighbors.resize(numNonConstNeighbors);
  auto index = 0;
  for (auto localj : mFlatToLocalIndex[locali]) {
    CHECK(size_t(localj) < mIsConstantBoundaryNode.size());
    if (!mIsConstantBoundaryNode[localj]) {
      CHECK(index < numNonConstNeighbors);
      CHECK(size_t(localj) < mLocalToGlobalIndex.size());
      globalNeighbors[index] = mLocalToGlobalIndex[localj];
      ++index;
    }
  }
  CHECK(index == numNonConstNeighbors);
}

template<typename Dimension>
void 
FlatConnectivity<Dimension>::
globalOverlapNeighborIndices(const int locali,
                             std::vector<int>& globalNeighbors) const {
  CHECK(mOverlapIndexingInitialized);
  CHECK(mBoundaryInformationInitialized);
  CHECK(mGlobalIndexingInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(size_t(locali) < mFlatOverlapToLocalIndex.size());

  const auto numNonConstNeighbors = mNumOverlapNeighbors[locali] - mNumConstantBoundaryOverlapNeighbors[locali];
  globalNeighbors.resize(numNonConstNeighbors);
  auto index = 0;
  for (auto localj : mFlatOverlapToLocalIndex[locali]) {
    CHECK(size_t(localj) < mIsConstantBoundaryNode.size());
    if (!mIsConstantBoundaryNode[localj]) {
      CHECK(index < numNonConstNeighbors);
      CHECK(size_t(localj) < mLocalToGlobalIndex.size());
      globalNeighbors[index] = mLocalToGlobalIndex[localj];
      ++index;
    }
  }
  CHECK(index == numNonConstNeighbors);
}

//------------------------------------------------------------------------------
// Get the local and global indices along with a map to make the global indices unique
//------------------------------------------------------------------------------
template<typename Dimension>
unsigned 
FlatConnectivity<Dimension>::
uniqueNeighborIndices(const unsigned locali,
                      std::vector<unsigned>& localNeighbors,
                      std::vector<unsigned>& globalNeighbors,
                      std::vector<unsigned>& constNeighbors,
                      std::vector<unsigned>& indexMap) const {
  TIME_FUNCTION;
  CHECK(mIndexingInitialized);
  CHECK(mBoundaryInformationInitialized);
  CHECK(mGlobalIndexingInitialized);
  CHECK(locali < static_cast<unsigned>(mNumConnectivityNodes));
  CHECK(size_t(locali) < mNumNeighbors.size() &&
        size_t(locali) < mFlatToLocalIndex.size());

  const auto numNonConstNeighbors = mNumNeighbors[locali] - mNumConstantBoundaryNeighbors[locali];
  const auto numConstNeighbors = mNumConstantBoundaryNeighbors[locali];

  // Get the original indices
  localNeighbors.resize(numNonConstNeighbors);
  globalNeighbors.resize(numNonConstNeighbors);
  constNeighbors.resize(numConstNeighbors);
  auto index = 0u;
  auto cindex = 0u;
  for (auto localj : mFlatToLocalIndex[locali]) {
    CHECK(size_t(localj) < mIsConstantBoundaryNode.size());
    if (mIsConstantBoundaryNode[localj]) {
      constNeighbors[cindex] = localj;
      ++cindex;
    }
    else {
      CHECK(index < static_cast<unsigned>(numNonConstNeighbors));
      CHECK(size_t(localj) < mLocalToGlobalIndex.size());
      localNeighbors[index] = localj;
      globalNeighbors[index] = mLocalToGlobalIndex[localj];
      ++index;
    }
  }
  CHECK(index == static_cast<unsigned>(numNonConstNeighbors));
  CHECK(cindex == static_cast<unsigned>(numConstNeighbors));

  // Make them unique
  return getUniqueIndices(globalNeighbors, indexMap);
}


template<typename Dimension>
unsigned
FlatConnectivity<Dimension>::
getUniqueIndices(const std::vector<unsigned>& globalNeighbors,
                 std::vector<unsigned>& indexMap) const {
  // if i is the index and ii is the unique index, indexMap[i] = ii
  const auto size = globalNeighbors.size();
  auto index = 0u;
  indexMap.resize(size);
  std::map<unsigned, unsigned> globalToIndex;
  for (auto j = 0u; j < size; ++j) {
    const auto global = globalNeighbors[j];
    auto it = globalToIndex.find(global);
    if (it == globalToIndex.end()) {
      // If the map doesn't have this global index yet, add a new index to the map
      globalToIndex[global] = index;
      indexMap[j] = index;
      ++index;
    }
    else {
      // If the map has this global index, then map back to that unique index
      indexMap[j] = it->second;
    }
  }

  return index;
}


//------------------------------------------------------------------------------
// Check whether NodeList ordering is appropriate for this function
//------------------------------------------------------------------------------
template<typename Dimension>
bool
FlatConnectivity<Dimension>::
fluidNodeListsFirst(const DataBase<Dimension>& dataBase) const
{
  auto nodeListi = 0;
  auto nodeListItri = dataBase.nodeListBegin();
  auto nodeListItrj = dataBase.fluidNodeListBegin();
  for (; nodeListItrj != dataBase.fluidNodeListEnd();
       ++nodeListItri, ++nodeListItrj, ++nodeListi) {
     if (*nodeListItri != *nodeListItrj)
     {
        return false;
     }
  }
  return true;
}

// //------------------------------------------------------------------------------
// // Check whether two points overlap, given the H values
// //------------------------------------------------------------------------------
// template<typename Dimension>
// bool
// FlatConnectivity<Dimension>::
// checkOverlap(const Vector& x1,
//              const SymTensor& H1,
//              const Vector& x2,
//              const SymTensor& H2,
//              const Scalar extent) const {
//   // For now, assume standard SPH
//   const auto h1 = Dimension::nDim / H1.Trace();
//   const auto h2 = Dimension::nDim / H2.Trace();

//   // Get the scaled distance
//   const auto x12Mag = (x1 - x2).magnitude();
//   const auto dist = x12Mag / (h1 + h2);

//   // Check whether scaled distance is less than kernel extent
//   return dist <= extent;
// }

} // end namespace Spheral
