//---------------------------------Spheral++----------------------------------//
// BilinearIndex
//
// Does a few simple integrals for testing
//----------------------------------------------------------------------------//
#include "BilinearIndex.hh"

#include "Geometry/CellFaceFlag.hh"
#include "DataBase/State.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
BilinearIndex<Dimension>::
BilinearIndex():
  mConnectivityInitialized(false),
  mConnectivityComputed(false),
  mSurfaceIndexingInitialized(false),
  mSurfaceIndexingComputed(false) {
}

//------------------------------------------------------------------------------
// Compute the connectivity information
//------------------------------------------------------------------------------
template<typename Dimension>
void
BilinearIndex<Dimension>::
computeConnectivity(const DataBase<Dimension>& dataBase) {
  if (!mConnectivityInitialized) {
    mBilinearFlatIndex = dataBase.newFluidFieldList(MapPair(), "bililnear flat index");
    mBilinearNodeIndex = dataBase.newFluidFieldList(VectorPair(), "bilinear node index");
    mConnectivityInitialized = true;
  }
  
  const auto numNodeLists = dataBase.numFluidNodeLists();
  const auto& connectivityMap = dataBase.connectivityMap();

  // Clear the index arrays and make a guess about how many indices there will be
  // Guess should work fine in 1d, but won't be as good in 2d/3d
  auto nodeListi = 0;
  for (auto nodeListItri = dataBase.fluidNodeListBegin();
       nodeListItri != dataBase.fluidNodeListEnd();
       ++nodeListItri, ++nodeListi) {
    const auto numNodesi = (*nodeListItri)->numNodes();
    for (auto nodei = 0u; nodei < numNodesi; ++nodei) {
      const auto numNeighbors = connectivityMap.numNeighborsForNode(nodeListi, nodei);
      mBilinearFlatIndex(nodeListi, nodei).clear();
      mBilinearFlatIndex(nodeListi, nodei).reserve(2 * numNeighbors);
      mBilinearNodeIndex(nodeListi, nodei).clear();
      mBilinearNodeIndex(nodeListi, nodei).reserve(2 * numNeighbors);
    }
  }
  
  // For each nodei, get the neighbors
  // Add the neighbors (nodej and nodek) to one another's connectivity
  // This forms an overlap connectivity, which may be the same as the connectivityMap one
  nodeListi = 0;
  for (auto nodeListItri = dataBase.fluidNodeListBegin();
       nodeListItri != dataBase.fluidNodeListEnd();
       ++nodeListItri, ++nodeListi) {
    const auto numNodesi = (*nodeListItri)->numInternalNodes(); // only over integration regions this process owns
    for (auto nodei = 0u; nodei < numNodesi; ++nodei) {
      const auto connectivity = connectivityMap.connectivityForNodeView(nodeListi, nodei);
      const auto pairi = std::make_pair(nodeListi, nodei);
      // Find overlap neighbors of the original point
      // For now, map the points to zero
      // We will index them once they are all added
      mBilinearFlatIndex(nodeListi, nodei).emplace(pairi, 0);
      for (auto nodeListj = 0u; nodeListj < numNodeLists; ++nodeListj) {
        for (auto nodej : connectivity[nodeListj]) {
          const auto pairj = std::make_pair(nodeListj, nodej);
          mBilinearFlatIndex(nodeListi, nodei).emplace(pairj, 0);
          mBilinearFlatIndex(nodeListj, nodej).emplace(pairi, 0);
          for (auto nodeListk = 0u; nodeListk < numNodeLists; ++nodeListk) {
            for (auto nodek : connectivity[nodeListk]) {
              const auto pairk = std::make_pair(nodeListk, nodek);
              mBilinearFlatIndex(nodeListj, nodej).emplace(pairk, 0);
            }
          }
        }
      }
    }
  }
  
  // Now that we have a unique set of point connectivity, we need to assign a unique index to each
  // We also create the connectivity information for the ghost nodes
  nodeListi = 0;
  for (auto nodeListItri = dataBase.fluidNodeListBegin();
       nodeListItri != dataBase.fluidNodeListEnd();
       ++nodeListItri, ++nodeListi) {
    const auto numNodesi = (*nodeListItri)->numNodes(); // includes ghost nodes
    for (auto nodei = 0u; nodei < numNodesi; ++nodei) {
      auto& flatIndex = mBilinearFlatIndex(nodeListi, nodei);
      const auto numIndices = flatIndex.size();
      auto& nodeIndex = mBilinearNodeIndex(nodeListi, nodei);
      nodeIndex.resize(numIndices);
      
      // For each set of flat indices, assign it a local index
      // and create the inverse of the index
      auto index = 0;
      for (auto it = flatIndex.begin(); it != flatIndex.end(); ++it, ++index) {
        it->second = index;
        nodeIndex[index] = it->first;
      }
    }
  }
  
  mConnectivityComputed = true;
}

//------------------------------------------------------------------------------
// Compute the surface connectivity information
//------------------------------------------------------------------------------
template<typename Dimension>
void
BilinearIndex<Dimension>::
computeSurfaceIndexing(const DataBase<Dimension>& dataBase,
                       const State<Dimension>& state) {
  if (!mSurfaceIndexingInitialized) {
    mSurfaceNormal = dataBase.newFluidFieldList(NormalType(), "surface normal");
    mSurfaceFlatIndex = dataBase.newFluidFieldList(MapNormal(), "local index for surface");
    mSurfaceFlags = dataBase.newFluidFieldList(std::vector<int>(), "surface flags");
    mSurfaceIndexingInitialized = true;
  }

  // Get DataBase values
  const auto numNodeLists = dataBase.numFluidNodeLists();
  const auto& connectivityMap = dataBase.connectivityMap();

  // Get State values
  const auto cells = state.fields(HydroFieldNames::cells, FacetedVolume());
  const auto cellFaceFlags = state.fields(HydroFieldNames::cellFaceFlags, std::vector<CellFaceFlag>());
  VERIFY(cells.size() == numNodeLists
         &&cellFaceFlags.size() == numNodeLists);
  
  // Clear the index arrays
  auto nodeListi = 0;
  for (auto nodeListItri = dataBase.fluidNodeListBegin();
       nodeListItri != dataBase.fluidNodeListEnd();
       ++nodeListItri, ++nodeListi) {
    const auto numNodesi = (*nodeListItri)->numNodes();
    for (auto nodei = 0u; nodei < numNodesi; ++nodei) {
      mSurfaceNormal(nodeListi, nodei).clear();
      mSurfaceFlatIndex(nodeListi, nodei).clear();
      mSurfaceFlags(nodeListi, nodei).clear();
    }
  }
  
  // Fill in the normals for each point
  ArrayDim normalArray;
  nodeListi = 0;
  for (auto nodeListItri = dataBase.fluidNodeListBegin();
       nodeListItri != dataBase.fluidNodeListEnd();
       ++nodeListItri, ++nodeListi) {
    const auto numNodesi = (*nodeListItri)->numNodes();
    for (auto nodei = 0u; nodei < numNodesi; ++nodei) {
      const auto& flags = cellFaceFlags(nodeListi, nodei);
      const auto numFlags = flags.size();
      if (numFlags > 0) {
        // Get the connectivity and surface information for this cell
        const auto connectivity = connectivityMap.connectivityForNodeView(nodeListi, nodei);
        const auto& cell = cells(nodeListi, nodei);
        const auto& facets = cell.facets();
        
        // If the surface flag is a void, then add the normal to the indexing
        for (auto flag : flags) {
          if (flag.nodeListj == -1) { // opposite side is a void
            mSurfaceFlags(nodeListi, nodei).push_back(flag.cellFace);
            const auto& facet = facets[flag.cellFace];
            const auto& normalArea = facet.normal();
            const auto normal = normalArea.unitVector();

            normalToArray(normal, normalArray);
            for (auto nodeListj = 0u; nodeListj < numNodeLists; ++nodeListj) {
              for (auto nodej : connectivity[nodeListj]) {
                // If the normal is added, the index is the current size of the surface normal
                const auto index = mSurfaceNormal(nodeListj, nodej).size();
                const auto added = mSurfaceFlatIndex(nodeListj, nodej).emplace(normalArray, index);
                if (added.second) {
                  // We avoid roundoff error by adding the normal directly instead of recreating it from the indices
                  mSurfaceNormal(nodeListj, nodej).push_back(normal);
                }
              }
            }

            // Add the self-contribution
            const auto index = mSurfaceNormal(nodeListi, nodei).size();
            const auto added = mSurfaceFlatIndex(nodeListi, nodei).emplace(normalArray, index);
            if (added.second) {
              mSurfaceNormal(nodeListi, nodei).push_back(normal);
            }
          }
        }
      }
    }
  }
  
  mSurfaceIndexingComputed = true;

  BEGIN_CONTRACT_SCOPE
  {
    nodeListi = 0;
    for (auto nodeListItri = dataBase.fluidNodeListBegin();
         nodeListItri != dataBase.fluidNodeListEnd();
         ++nodeListItri, ++nodeListi) {
      const auto numNodesi = (*nodeListItri)->numNodes();
      for (auto nodei = 0u; nodei < numNodesi; ++nodei) {
        const auto ni = std::make_pair(nodeListi, nodei);

        // Make sure the sizes of the map and normals match up
        CHECK(mSurfaceNormal(nodeListi, nodei).size() == mSurfaceFlatIndex(nodeListi, nodei).size());
        
        // Make sure the map returns the expected surface index
        const auto numSurfaces = numSurfacesForRow(ni);
        const auto& norm = normal(ni);
        for (auto s = 0; s < numSurfaces; ++s) {
          normalToArray(norm[s], normalArray);
          CHECK(s == surfaceIndex(ni, normalArray));
        }

        // Make sure all key/value pairs line up to an existing normal
        for (auto& mapVal : mSurfaceFlatIndex(nodeListi, nodei)) {
          const auto& normKey = mapVal.first;
          CONTRACT_VAR(normKey);
          const auto s = mapVal.second;
          CHECK(s < numSurfaces);
          normalToArray(norm[s],
                        normalArray);
          CHECK(normalArray == normKey);
        }
      }
    }
  }
  END_CONTRACT_SCOPE
}

} // end namespace Spheral
