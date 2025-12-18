//------------------------------------------------------------------------------
// Generate a mesh for the given set of NodeLists.
//------------------------------------------------------------------------------
#include "generateMesh.hh"
#include "computeGenerators.hh"
#include "Mesh.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/generateVoidNodes.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/Timer.hh"

#include <algorithm>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

template<typename Dimension, typename NodeListIterator, typename BoundaryIterator>
void
generateMesh(const NodeListIterator nodeListBegin,
             const NodeListIterator nodeListEnd,
             const BoundaryIterator boundaryBegin,
             const BoundaryIterator boundaryEnd,
             const typename Dimension::Vector& xmin,
             const typename Dimension::Vector& xmax,
             const bool meshGhostNodes,
             const bool generateVoid,
             const bool /*generateParallelConnectivity*/,
             const bool removeBoundaryZones,
             const double voidThreshold,
             Mesh<Dimension>& mesh,
             NodeList<Dimension>& voidNodes) {

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  // The total number of NodeLists we're working on.
  const size_t numNodeLists = distance(nodeListBegin, nodeListEnd);
  CONTRACT_VAR(numNodeLists);
  const unsigned voidOffset = distance(nodeListBegin, find(nodeListBegin, nodeListEnd, &voidNodes));

  // Pre-conditions.
  VERIFY2(voidOffset == numNodeLists - 1,
          "You must ensure the void nodes are part of the node set (at the end).");
  VERIFY2(voidNodes.numInternalNodes() == 0,
          "Void nodes not empty on startup!");
  VERIFY2(not (generateVoid and removeBoundaryZones),
          "You cannot simultaneously request generateVoid and removeBoundaryZones.");

  // Extract the set of generators this domain needs.
  // This method gives us both the positions and Hs for the generators.
  TIME_BEGIN("generateMesh::computeGenerators");
  vector<Vector> generators;
  vector<SymTensor> Hs;
  vector<unsigned> offsets;
  computeGenerators<Dimension, NodeListIterator, BoundaryIterator>(nodeListBegin, nodeListEnd, 
                                                                   boundaryBegin, boundaryEnd,
                                                                   meshGhostNodes,
                                                                   xmin, xmax, 
                                                                   generators, Hs, offsets);
  TIME_END("generateMesh::computeGenerators");

  // Construct the mesh.
  TIME_BEGIN("generateMesh::reconstruct");
  mesh.reconstruct(generators, xmin, xmax, boundaryBegin, boundaryEnd);
  CHECK(mesh.numZones() == generators.size());
  TIME_END("generateMesh::reconstruct");

  // Are we generating void?
  if (generateVoid or removeBoundaryZones) {
    TIME_BEGIN("generateMesh::generateVoidNodes");
    unsigned numInternal = 0;
    double nPerh = 0;
    for (NodeListIterator itr = nodeListBegin; itr != nodeListEnd - 1; ++itr) {
      numInternal += (**itr).numInternalNodes();
      nPerh = (**itr).nodesPerSmoothingScale();
    }
    mesh.generateParallelRind(generators, Hs);
    generateVoidNodes(generators, Hs, mesh, xmin, xmax, numInternal, nPerh, voidThreshold, voidNodes);

    computeGenerators<Dimension, NodeListIterator, BoundaryIterator>(nodeListBegin, nodeListEnd, 
                                                                     boundaryBegin, boundaryEnd,
                                                                     meshGhostNodes,
                                                                     xmin, xmax, 
                                                                     generators, Hs, offsets);
    TIME_END("generateMesh::generateVoidNodes");
    // Construct the mesh.
    TIME_BEGIN("generateMesh::generateVoidNodes::reconstruct");
    mesh.reconstruct(generators, xmin, xmax, boundaryBegin, boundaryEnd);
    TIME_END("generateMesh::generateVoidNodes::reconstruct");
    CHECK(mesh.numZones() == generators.size());

  }

  // Remove any zones for generators that are not local to this domain.
  TIME_BEGIN("generateMesh::cullZones");
  vector<unsigned> mask(mesh.numZones(), 0);
  unsigned ioff = 0;
  for (NodeListIterator itr = nodeListBegin; itr != nodeListEnd; ++itr, ++ioff) {
    if (ioff != voidOffset) {
      const unsigned zoneOffset = offsets[ioff];
      if (meshGhostNodes) {
        fill(mask.begin() + zoneOffset, mask.begin() + offsets[ioff + 1], 1);
      } else {
        fill(mask.begin() + zoneOffset, mask.begin() + zoneOffset + (*itr)->numInternalNodes(), 1);
      }
    }
  }
  if (not removeBoundaryZones) {
    fill(mask.begin() + offsets[voidOffset], mask.begin() + offsets[voidOffset] + voidNodes.numInternalNodes(), 1);
  }
  mesh.removeZonesByMask(mask);

  // If we removed boundary nodes, then the void was created for this purpose
  // alone.
  if (removeBoundaryZones) {
    voidNodes.numInternalNodes(0);
    offsets.back() = mesh.numZones();
  }
  TIME_END("generateMesh::cullZones");

  // Fill in the offset information.
  mesh.storeNodeListOffsets(nodeListBegin, nodeListEnd, offsets);

  // That's it.
}

}
