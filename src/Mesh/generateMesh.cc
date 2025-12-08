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
#include "Utilities/timingUtilities.hh"

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
  // if (Process::getRank() == 0) cerr << "Computing generators" << endl;
  // Timing::Time t0 = Timing::currentTime();
  vector<Vector> generators;
  vector<SymTensor> Hs;
  vector<unsigned> offsets;
  computeGenerators<Dimension, NodeListIterator, BoundaryIterator>(nodeListBegin, nodeListEnd, 
                                                                   boundaryBegin, boundaryEnd,
                                                                   meshGhostNodes,
                                                                   xmin, xmax, 
                                                                   generators, Hs, offsets);
  // if (Process::getRank() == 0) cerr << "generateMesh:: required " 
  //                                   << Timing::difference(t0, Timing::currentTime())
  //                                   << " seconds to construct generators." << endl;

  // Construct the mesh.
  // t0 = Timing::currentTime();
  mesh.reconstruct(generators, xmin, xmax, boundaryBegin, boundaryEnd);
  CHECK(mesh.numZones() == generators.size());
  // if (Process::getRank() == 0) cerr << "generateMesh:: required " 
  //                                   << Timing::difference(t0, Timing::currentTime())
  //                                   << " seconds to construct mesh." << endl;

  // Are we generating void?
  // t0 = Timing::currentTime();
  if (generateVoid or removeBoundaryZones) {
    // if (Process::getRank() == 0)  cerr << "Computing void nodes." << endl;
    unsigned numInternal = 0;
    double nPerh = 0;
    for (NodeListIterator itr = nodeListBegin; itr != nodeListEnd - 1; ++itr) {
      numInternal += (**itr).numInternalNodes();
      nPerh = (**itr).nodesPerSmoothingScale();
    }
    mesh.generateParallelRind(generators, Hs);
    generateVoidNodes(generators, Hs, mesh, xmin, xmax, numInternal, nPerh, voidThreshold, voidNodes);

    // if (Process::getRank() == 0) cerr << "Recomputing generators with void." << endl;
    computeGenerators<Dimension, NodeListIterator, BoundaryIterator>(nodeListBegin, nodeListEnd, 
                                                                     boundaryBegin, boundaryEnd,
                                                                     meshGhostNodes,
                                                                     xmin, xmax, 
                                                                     generators, Hs, offsets);
    // if (Process::getRank() == 0) cerr << "generateMesh:: required " 
    //                                   << Timing::difference(t0, Timing::currentTime())
    //                                   << " seconds to construct generators." << endl;

    // Construct the mesh.
    // t0 = Timing::currentTime();
    mesh.reconstruct(generators, xmin, xmax, boundaryBegin, boundaryEnd);
    CHECK(mesh.numZones() == generators.size());
    // if (Process::getRank() == 0) cerr << "generateMesh:: required " 
    //                                   << Timing::difference(t0, Timing::currentTime())
    //                                   << " seconds to construct mesh." << endl;

  }

  // Remove any zones for generators that are not local to this domain.
  // if (Process::getRank() == 0) cerr << "Removing zones" << endl;
  // t0 = Timing::currentTime();
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
  // if (Process::getRank() == 0) cerr << "generateMesh:: required " 
  //                                   << Timing::difference(t0, Timing::currentTime())
  //                                   << " seconds to remove boundary elements." << endl;

  // // If requested we also compute the parallel connectivity.
  // // if (Process::getRank() == 0) cerr << "Computing parallel connectivity" << endl;
  // if (generateParallelConnectivity) {
  //   t0 = Timing::currentTime();
  //   mesh.generateDomainInfo();
  //   if (Process::getRank() == 0) cerr << "generateMesh:: required " 
  //                                     << Timing::difference(t0, Timing::currentTime())
  //                                     << " seconds to generate parallel connectivity." << endl;
  // }

  // Fill in the offset information.
  mesh.storeNodeListOffsets(nodeListBegin, nodeListEnd, offsets);

  // That's it.
}

}
