//------------------------------------------------------------------------------
// Look for any points that touch a surface (multi-material or void).
// For such points:
//   - Remove any non-surface multi-material overlap.
//   - If not a surface point, flag this point as touching a surface point with
//     surfacePoint=-1.
//------------------------------------------------------------------------------
#include "editMultimaterialSurfaceTopology.hh"
#include "Utilities/DBC.hh"
#include "Utilities/Timer.hh"

#include <vector>
#include <algorithm>

using std::vector;
using std::cerr;
using std::endl;

namespace Spheral {

template<typename Dimension>
void
editMultimaterialSurfaceTopology(FieldList<Dimension, int>& surfacePoint,
                                 ConnectivityMap<Dimension>& connectivityMap) {
  TIME_BEGIN("CRKSPH_editMultimaterialSurfaceTopology");

  // Declare some useful stuff and preconditions.
  const auto  numNodeLists = surfacePoint.size();
  const auto& nodeLists = connectivityMap.nodeLists();
  REQUIRE(nodeLists.size() == numNodeLists);

  // Record the points we want to remove from the neighbor set.
  FieldList<Dimension, vector<vector<int>>> neighborsToCut(FieldStorageType::CopyFields);
  for (auto nodeList: nodeLists) neighborsToCut.appendNewField("cut neighbors", *nodeList, vector<vector<int>>(numNodeLists));

  // Here we go.
  for (auto iNodeList = 0u; iNodeList < numNodeLists - 1; ++iNodeList) {
    const auto n = nodeLists[iNodeList]->numInternalNodes();

#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto allneighbors = connectivityMap.connectivityForNodeView(iNodeList, i);

      // printf(" --> (%d, %d) :", iNodeList, i);

      if (surfacePoint(iNodeList, i) == 0) {

        // This is not a surface point: cut all connections to other NodeLists.
        for (auto jNodeList = 0u; jNodeList < numNodeLists; ++jNodeList) {
          if (jNodeList != iNodeList) {
            neighborsToCut(iNodeList, i)[jNodeList].resize(allneighbors[jNodeList].size());
            for (auto k = 0u; k < allneighbors[jNodeList].size(); ++k) neighborsToCut(iNodeList, i)[jNodeList][k] = k;
          }
        }

      } else {

        // This is a surface point, so we allow connectivity to surface points in other NodeLists.
        for (auto jNodeList = 0u; jNodeList < numNodeLists; ++jNodeList) {
          if (jNodeList != iNodeList) {
            const bool ijboundary = (surfacePoint(iNodeList, i) && (1 << (jNodeList + 1)) > 0);
            if (ijboundary) {
              // i is a neighbor to jNodeList, so we need to pick through the points and find like neighbors.
              const auto& neighbors = allneighbors[jNodeList];
              for (auto k = 0u; k < neighbors.size(); ++k) {
                const auto j = neighbors[k];
                const bool jiboundary = (surfacePoint(jNodeList, j) && (1 << (iNodeList + 1)) > 0);
                if (not jiboundary) neighborsToCut(iNodeList, i)[jNodeList].push_back(k);
              }
            } else {
              // Since i is not a boundary point to this neighbor NodeList, cut all ties.
              neighborsToCut(iNodeList, i)[jNodeList].resize(allneighbors[jNodeList].size());
              for (auto k = 0u; k < allneighbors[jNodeList].size(); ++k) neighborsToCut(iNodeList, i)[jNodeList][k] = k;
            }
          }
        }
      }
    }
  }

  // Apply our toplogical cuts to the ConnectivityMap.
  connectivityMap.removeConnectivity(neighborsToCut);

  // // Now check for any non-surface nodes that overlap surface in the newly reduced topology.
  // // Flag those specially with -1 in surfacePoint.
  // for (auto iNodeList = 0; iNodeList < numNodeLists; ++iNodeList) {
  //   const auto n = nodeLists[iNodeList]->numInternalNodes();
  //   for (auto i = 0; i < n; ++i) {
  //     if (surfacePoint(iNodeList, i) == 0) {
  //       const auto& allneighbors = connectivityMap.connectivityForNode(iNodeList, i);
  //       const auto& neighbors = allneighbors[iNodeList];
  //       for (auto j: neighbors) {
  //         if (surfacePoint(iNodeList, j) > 0) {
  //           surfacePoint(iNodeList, i) = -1;
  //           break;
  //         }
  //       }
  //     }
  //   }
  // }

  // // Go out one more rind of points two-steps removed from the boundary.
  // for (auto iNodeList = 0; iNodeList < numNodeLists; ++iNodeList) {
  //   const auto n = nodeLists[iNodeList]->numInternalNodes();
  //   for (auto i = 0; i < n; ++i) {
  //     if (surfacePoint(iNodeList, i) == -2) {
  //       const auto& allneighbors = connectivityMap.connectivityForNode(iNodeList, i);
  //       const auto& neighbors = allneighbors[iNodeList];
  // 	for (auto j: neighbors) {
  // 	  if (surfacePoint(iNodeList, j) == 0) surfacePoint(iNodeList, j) = -3;
  // 	}
  //     }
  //   }
  // }

  // // Flip any -2 points to -1.
  // for (auto fieldPtr: surfacePoint) {
  //   const auto n = fieldPtr->numInternalElements();
  //   std::replace(fieldPtr->begin(), fieldPtr->begin() + n, -2, -1);
  // }

  // // BLAGO!
  // for (auto fieldPtr: surfacePoint) {
  //   const auto n = fieldPtr->numInternalElements();
  //   const auto& pos = fieldPtr->nodeListPtr()->positions();
  //   for (auto i = 0; i < n; ++i) {
  //     if ((*fieldPtr)(i) == 0 and std::abs(pos(i).x()) < 0.25) (*fieldPtr)(i) = -3;
  //   }
  // }
  // // BLAGO!

  TIME_END("CRKSPH_editMultimaterialSurfaceTopology");
}

}
