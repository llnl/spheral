//---------------------------------Spheral++----------------------------------//
// PolyhedralMesh -- 3-D mesh class.
//
// Created by JMO, Wed Jan  5 21:09:11 PST 2011
//----------------------------------------------------------------------------//
#ifdef USE_POLYTOPE
#include "polytope/polytope.hh"
#endif

#include "Mesh.hh"
#include "MeshConstructionUtilities.hh"
#include "Utilities/DBC.hh"
#include "Utilities/timingUtilities.hh"
#include "Distributed/Communicator.hh"

#include <limits>
#include <set>
#include <sstream>
using std::vector;
using std::map;
using std::set;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {


//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<> const unsigned Mesh<Dim<3> >::minFacesPerZone = 4;
template<> const unsigned Mesh<Dim<3> >::minEdgesPerZone = 6;
template<> const unsigned Mesh<Dim<3> >::minNodesPerZone = 4;
template<> const unsigned Mesh<Dim<3> >::minEdgesPerFace = 3;
template<> const unsigned Mesh<Dim<3> >::minNodesPerFace = 3;

//------------------------------------------------------------------------------
// Mesh::reconstructInternal
//------------------------------------------------------------------------------
template<>
void
Mesh<Dim<3> >::
reconstructInternal(const vector<Dim<3>::Vector>& generators,
                    const Dim<3>::Vector& xmin,
                    const Dim<3>::Vector& xmax) {
  CONTRACT_VAR(generators);
  CONTRACT_VAR(xmin);
  CONTRACT_VAR(xmax);
#ifdef USE_POLYTOPE

  // Some useful typedefs.
  typedef Dim<3> Dimension;
  typedef pair<unsigned, unsigned> EdgeHash;

  // Is there anything to do?
  if (generators.size() == 0) return;

  // The tolerance on which we consider positions to be degenerate.
  const double xtol = 1.0e-12*(xmax - xmin).maxElement();
  const double xtol2 = xtol*xtol;

  // Pre-conditions.
  int i, j, k, igen, jgen;
  const unsigned numGens = generators.size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(xmin.x() < xmax.x() and
            xmin.y() < xmax.y() and
            xmin.z() < xmax.z());
    for (igen = 0; igen != numGens; ++igen)
      REQUIRE2(testPointInBox(generators[igen], xmin, xmax), 
               "Generator out of bounds:  " << generators[igen] << " not in [" << xmin << " " << xmax << "]");
    for (i = 0; i < numGens - 1; ++i) {
      for (j = i + 1; j < numGens; ++j) {
        REQUIRE2((generators[i] - generators[j]).magnitude2() > xtol2, 
                 "Degenerate generator positions:  " << i << " " << j << " " << generators[i] << " " << generators[j]);
      }
    }
  }
  END_CONTRACT_SCOPE

  // The inverse box scale.
  const Vector box = xmax - xmin;
  const Vector boxInv(safeInv(box.x()),
                      safeInv(box.y()),
                      safeInv(box.z()));

  // Build the polytope style generator positions as a single flat vector.
  vector<double> gens;
  gens.reserve(3*generators.size());
  for (igen = 0; igen != numGens; ++igen) {
    gens.push_back(generators[igen].x());
    gens.push_back(generators[igen].y());
    gens.push_back(generators[igen].z());
  }
  CHECK(gens.size() == 3*numGens);

  // Do the polytope tessellation.  We use the TetGen based tessellator.
  Timing::Time t0 = Timing::currentTime();
  polytope::Tessellation<3, double> tessellation;
//   {
// #ifdef SPHERAL_ENABLE_MPI
//     polytope::SerialDistributedTessellator<3, double> tessellator
// #if defined USE_TETGEN && ( USE_TETGEN>0 )
//         (new polytope::TetgenTessellator(),
// #else
//         (new polytope::VoroPP_3d<double>(),
// #endif
//          true,     // Manage memory for serial tessellator
//          true);    // Build parallel connectivity
// #else  // not ENABLE_MPI
// #if defined USE_TETGEN && ( USE_TETGEN>0 )
//     polytope::TetgenTessellator tessellator ;
// #else
//     polytope::VoroPP_3d<double> tessellator ;
// #endif
// #endif  // ENABLE_MPI

//     // Bounded Voronoi tessellation
//     tessellator.tessellate(gens, 
//                            const_cast<double*>(xmin.begin()), 
//                            const_cast<double*>(xmax.begin()),
//                            tessellation);
// }
  CHECK(tessellation.cells.size() == numGens);
  if (Process::getRank() == 0) cerr << "PolyhedralMesh:: required " 
                                    << Timing::difference(t0, Timing::currentTime())
                                    << " seconds to construct polytope tessellation." << endl;

  // Read out and construct the Mesh node positions.
  t0 = Timing::currentTime();
  const unsigned numNodes = tessellation.nodes.size()/3;
  mNodePositions.reserve(numNodes);
  for (i = 0; i != numNodes; ++i) {
    mNodePositions.push_back(Vector(max(xmin.x(), min(xmax.x(), xmin.x() + tessellation.nodes[3*i]  )),
                                    max(xmin.y(), min(xmax.y(), xmin.y() + tessellation.nodes[3*i+1])),
                                    max(xmin.z(), min(xmax.z(), xmin.z() + tessellation.nodes[3*i+2]))));
  }

  // Build the edges and faces.
  const unsigned numFaces = tessellation.faces.size();
  unsigned inode, jnode, iedge, n, numEdges = 0;
  EdgeHash ehash;
  map<EdgeHash, unsigned> edgeHash2ID;
  map<EdgeHash, unsigned>::iterator emapItr;
  map<unsigned, set<unsigned> > nodeZones;
  for (i = 0; i != numFaces; ++i) {
    n = tessellation.faces[i].size();
    CHECK(n >= 3);
    CHECK(tessellation.faceCells[i].size() == 1 or
          tessellation.faceCells[i].size() == 2);
    igen = tessellation.faceCells[i][0];
    jgen = (tessellation.faceCells[i].size() == 2 ? tessellation.faceCells[i][1] :
            igen < 0 ? UNSETID : ~UNSETID);
    vector<unsigned> faceEdges;
    for (j = 0; j != n; ++j) {
      inode = tessellation.faces[i][j];
      jnode = tessellation.faces[i][(j + 1) % n];
      ehash = hashEdge(inode, jnode);
      emapItr = edgeHash2ID.find(ehash);
      if (emapItr == edgeHash2ID.end()) {
        iedge = numEdges;
        edgeHash2ID[ehash] = numEdges++;
        mEdges.push_back(Edge(*this, iedge, inode, jnode));
      } else {
        iedge = emapItr->second;
      }
      faceEdges.push_back(iedge);
    }
    mFaces.push_back(Face(*this, i, igen, jgen, faceEdges));
    for (int j: tessellation.faceCells[i]) {
      j = positiveID(j);
      nodeZones[inode].insert(j);
      nodeZones[jnode].insert(j);
    }
  }
  CHECK(mEdges.size() == numEdges);
  CHECK(mFaces.size() == numFaces);

  // Construct the nodes.
  mNodes.reserve(numNodes);
  for (i = 0; i != numNodes; ++i) mNodes.push_back(Node(*this, i, vector<unsigned>(nodeZones[i].begin(),
                                                                                   nodeZones[i].end())));
  CHECK(mNodes.size() == numNodes);

  // Construct the zones.
  for (i = 0; i != numGens; ++i) mZones.push_back(Zone(*this, i, tessellation.cells[i]));
  CHECK(mZones.size() == numGens);

  // Copy the parallel info.
  mNeighborDomains = tessellation.neighborDomains;
  mSharedNodes = tessellation.sharedNodes;
  mSharedFaces = tessellation.sharedFaces;

  // Report our final timing and we're done.
  if (Process::getRank() == 0) cerr << "PolyhedralMesh:: required " 
                                    << Timing::difference(t0, Timing::currentTime())
                                    << " seconds to construct mesh elements." << endl;

#endif // USE_POLYTOPE
}

//------------------------------------------------------------------------------
// Mesh::reconstructInternal (FacetedVolume)
//------------------------------------------------------------------------------
template<>
void
Mesh<Dim<3> >::
reconstructInternal(const vector<Dim<3>::Vector>& /*generators*/,
                    const Dim<3>::FacetedVolume& /*boundary*/) {
  VERIFY2(false, "PolyhedralMesh ERROR: tessellations with arbitrary polyhedral boundaries not supported.");
}

//------------------------------------------------------------------------------
// Compute the bounding surface of the mesh.
//------------------------------------------------------------------------------
template<>
Dim<3>::FacetedVolume
Mesh<Dim<3> >::
boundingSurface() const {

  // Flatten the set of communicated nodes into a set.
  set<unsigned> sharedNodes;
  unsigned domainID;
  for (domainID = 0; domainID != mSharedNodes.size(); ++domainID) {
    std::copy(mSharedNodes[domainID].begin(), mSharedNodes[domainID].end(),
              std::inserter(sharedNodes, sharedNodes.begin()));
  }

  // Build the global IDs for the mesh nodes.
  const vector<unsigned> local2globalIDs = this->globalMeshNodeIDs();

  // Look for the faces that bound the mesh.  We build up the global
  // vertex indices, and the associated positions.
  bool useFace;
  unsigned i, j, iglobal, izone;
  map<unsigned, Vector> globalVertexPositions;
  vector<vector<unsigned> > facetIndices;
  for (const Face& face: mFaces) {
    const vector<unsigned>& nodeIDs = face.nodeIDs();
    useFace = (face.zone1ID() == (int)UNSETID or face.zone2ID() == (int)UNSETID);
    i = 0;
    while (useFace and i != nodeIDs.size()) {
      useFace = (sharedNodes.find(nodeIDs[i]) == sharedNodes.end());
      ++i;
    }
    if (useFace) {
      vector<unsigned> ids;
      for (const unsigned i: nodeIDs) {
        CHECK(i < local2globalIDs.size());
        iglobal = local2globalIDs[i];
        CHECK(globalVertexPositions.find(iglobal) == globalVertexPositions.end() or
              fuzzyEqual((globalVertexPositions[iglobal] - mNodePositions[i]).magnitude2(), 0.0));
        globalVertexPositions[iglobal] = mNodePositions[i];
        ids.push_back(iglobal);
      }

      // Do we need to reverse the face orientation?
      izone = std::min(face.zone1ID(), face.zone2ID());
      CHECK(izone != UNSETID);
      if (face.compare(mZones[izone].position()) == 1) std::reverse(ids.begin(), ids.end());

      facetIndices.push_back(ids);
    }
  }

#ifdef SPHERAL_ENABLE_MPI
  // In the parallel case we have to construct the total surface and distribute
  // it to everyone.
  //const unsigned rank = Process::getRank();
  const unsigned numDomains = Process::getTotalNumberOfProcesses();
  if (numDomains > 1) {
    unsigned bufSize, nfacets;
    vector<char> localBuffer, buffer;
    vector<char>::const_iterator bufItr;
    vector<Vector> otherVertices;
    vector<unsigned> otherIndices;

    // Pack our local data, and then erase our local copy to be rebuilt 
    // consistently for everyone.
    packElement(globalVertexPositions, localBuffer);
    packElement(unsigned(facetIndices.size()), localBuffer);
    for (i = 0; i != facetIndices.size(); ++i) packElement(facetIndices[i], localBuffer);
    globalVertexPositions = map<unsigned, Vector>();
    facetIndices = vector<vector<unsigned> >();

    // Distribute the complete data to everyone.
    for (domainID = 0; domainID != numDomains; ++domainID) {
      buffer = localBuffer;
      bufSize = localBuffer.size();
      MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, domainID, Communicator::communicator());
      if (bufSize > 0) {
        buffer.resize(bufSize);
        MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, domainID, Communicator::communicator());
        bufItr = buffer.begin();
        unpackElement(globalVertexPositions, bufItr, buffer.end());
        unpackElement(nfacets, bufItr, buffer.end());
        for (i = 0; i != nfacets; ++i) {
          otherIndices = vector<unsigned>();
          unpackElement(otherIndices, bufItr, buffer.end());
          CHECK2(otherIndices.size() >= 3, "Bad size: " << otherIndices.size());

          facetIndices.push_back(otherIndices);
        }
      }
    }
  }
#endif

  // Extract the vertex positions as an array, and map the global IDs 
  // to index in this array.
  map<unsigned, unsigned> global2vertexID;
  vector<Vector> vertices;
  vertices.reserve(globalVertexPositions.size());
  i = 0;
  for (map<unsigned, Vector>::const_iterator itr = globalVertexPositions.begin();
       itr != globalVertexPositions.end();
       ++itr, ++i) {
    global2vertexID[itr->first] = i;
    vertices.push_back(itr->second);
  }
  CHECK(i == globalVertexPositions.size());
  CHECK(global2vertexID.size() == globalVertexPositions.size());
  CHECK(vertices.size() == globalVertexPositions.size());

  // Transform the facet node indices to the vertex array numbering.
  for (vector<unsigned>& indices: facetIndices) {
    CHECK(indices.size() >= 3);
    for (j = 0; j != indices.size(); ++j) {
      CHECK(global2vertexID.find(indices[j]) != global2vertexID.end());
      indices[j] = global2vertexID[indices[j]];
      CHECK(indices[j] < vertices.size());
    }
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    for (const vector<unsigned>& indices: facetIndices) {
      CONTRACT_VAR(indices);
      ENSURE(indices.size() >= 3);
      ENSURE(*max_element(indices.begin(), indices.end()) < vertices.size());
    }
  }
  END_CONTRACT_SCOPE

  // That's it.
  return FacetedVolume(vertices, facetIndices);
}

//------------------------------------------------------------------------------
// Internal method add new mesh elements for existing node positions.
//------------------------------------------------------------------------------
template<>
void
Mesh<Dim<3> >::
createNewMeshElements(const vector<vector<vector<unsigned> > >& newCells) {

  typedef pair<unsigned, unsigned> EdgeHash;
  typedef set<unsigned> FaceHash;
  typedef pair<int, int> FaceZoneHash;

  // Pre-conditions.
  REQUIRE(mNodes.size() <= mNodePositions.size());
  BEGIN_CONTRACT_SCOPE
  {
    for (const vector<vector<unsigned> >& cellFaces: newCells) {
      REQUIRE(cellFaces.size() >= minFacesPerZone);
      for (const vector<unsigned>& faceNodes: cellFaces) {
        REQUIRE(faceNodes.size() >= minNodesPerFace);
        for (const unsigned inode: faceNodes) {
          CONTRACT_VAR(inode);
          REQUIRE(inode < mNodePositions.size());
        }
      }
    }
  }
  END_CONTRACT_SCOPE

  // Some useful sizes.
  const unsigned numOldNodes = mNodes.size();
  const unsigned numNewNodes = mNodePositions.size();
  const unsigned numOldEdges = mEdges.size();
  const unsigned numOldFaces = mFaces.size();
  const unsigned numOldZones = mZones.size();

  // Copy the existing face->zone connectivity.
  map<FaceHash, FaceZoneHash> faceZones;
  for (unsigned iface = 0; iface != numOldFaces; ++iface) {
    faceZones[FaceHash(mFaces[iface].mNodeIDs.begin(), mFaces[iface].mNodeIDs.end())] = FaceZoneHash(mFaces[iface].mZone1ID, mFaces[iface].mZone2ID);
  }

  // Add any new face->zone elements.
  for (unsigned i = 0; i != newCells.size(); ++i) {
    const vector<vector<unsigned> >& zoneFaces = newCells[i];
    const int newZoneID = numOldZones + i;
    // cerr << "New zone " << newZoneID << " with nodes : ";
    // for (const vector<unsigned>& faceNodes: zoneFaces) {
    //   copy(faceNodes.begin(), faceNodes.end(), ostream_iterator<unsigned>(cerr, " "));
    //   cerr << " : ";
    // }
    // cerr << endl;
    for (const vector<unsigned>& faceNodes: zoneFaces) {
      // cerr << "Adding face with nodes : ";
      // copy(faceNodes.begin(), faceNodes.end(), ostream_iterator<unsigned>(cerr, " "));
      const FaceHash fhash(faceNodes.begin(), faceNodes.end());
      map<FaceHash, FaceZoneHash>::iterator faceItr = faceZones.find(fhash);
      if (faceItr == faceZones.end()) {
        // New face.
        faceZones[fhash] = FaceZoneHash(newZoneID, ~UNSETID);
      } else {
        // Existing face, which means one of the cells for this face had better be UNSETID.
        FaceZoneHash& zones = faceItr->second;
        // if (!(positiveID(zones.first) == UNSETID or positiveID(zones.second) == UNSETID)) {
        //   cerr << positiveID(zones.first) << " : ";
        //   copy(
        CHECK2(positiveID(zones.first) == UNSETID or positiveID(zones.second) == UNSETID,
               zones.first << " " << zones.second << " : " << numOldZones << " " << newZoneID);
        if (positiveID(zones.first) == UNSETID) {
          zones.first = zones.first < 0 ? ~newZoneID : newZoneID;
        } else {
          zones.second = zones.second < 0 ? ~newZoneID : newZoneID;
        }
      }
      // cerr << "  :  " << faceZones[fhash].first << " " << faceZones[fhash].second << endl;
    }
  }

  // Based on the face->zone connectivity we reconstruct the node->zone connectivity.
  map<unsigned, set<unsigned> > nodeZones;
  for (map<FaceHash, FaceZoneHash>::const_iterator faceItr = faceZones.begin();
       faceItr != faceZones.end();
       ++faceItr) {
    const FaceHash& nodes = faceItr->first;
    const unsigned zone1 = positiveID(faceItr->second.first);
    const unsigned zone2 = positiveID(faceItr->second.second);
    for (const unsigned inode: nodes) {
      nodeZones[inode].insert(zone1);
      nodeZones[inode].insert(zone2);
    }
  }

  // Update the node->zones for existing nodes.
  for (unsigned inode = 0; inode != numOldNodes; ++inode) {
    mNodes[inode].mZoneIDs = vector<unsigned>(nodeZones[inode].begin(), nodeZones[inode].end());
  }

  // Create the new nodes.
  mNodes.reserve(numNewNodes);
  for (unsigned inode = numOldNodes; inode != numNewNodes; ++inode) {
    mNodes.push_back(Node(*this, inode, vector<unsigned>(nodeZones[inode].begin(), nodeZones[inode].end())));
  }
  CHECK(mNodes.size() == numNewNodes);

  // Determine the existing edge hash->edgeID mapping.
  map<EdgeHash, unsigned> edgeHash2ID;
  for (unsigned iedge = 0; iedge != numOldEdges; ++iedge) {
    edgeHash2ID[hashEdge(mEdges[iedge].mNode1ID, mEdges[iedge].mNode2ID)] = iedge;
  }

  // Similarly get the existing face hash->faceID mapping, hashing based on the face nodes.
  // We simultaneously update the face->zone connectivity.
  map<FaceHash, unsigned> faceHash2ID;
  for (unsigned iface = 0; iface != numOldFaces; ++iface) {
    const FaceHash fhash(mFaces[iface].mNodeIDs.begin(), mFaces[iface].mNodeIDs.end());
    const FaceZoneHash& zones = faceZones[fhash];
    faceHash2ID[fhash] = iface;
    mFaces[iface].mZone1ID = zones.first;
    mFaces[iface].mZone2ID = zones.second;
  }

  // Create any new edges, faces, and zones.
  for (unsigned i = 0; i != newCells.size(); ++i) {
    const int newZoneID = numOldZones + i;
    const vector<vector<unsigned> >& zoneFaceNodes = newCells[i];
    vector<int> zoneFaces;
    for (const vector<unsigned>& faceNodes: zoneFaceNodes) {
      vector<unsigned> faceEdges;
      const unsigned n = faceNodes.size();
      for (unsigned k = 0; k != n; ++k) {
        const unsigned inode1 = faceNodes[k];
        const unsigned inode2 = faceNodes[(k + 1) % n];
        const EdgeHash ehash = hashEdge(inode1, inode2);
        int iedge;
        const map<EdgeHash, unsigned>::const_iterator itr = edgeHash2ID.find(ehash);
        if (itr == edgeHash2ID.end()) {
          iedge = edgeHash2ID.size();
          edgeHash2ID[ehash] = iedge;
          mEdges.push_back(Edge(*this, iedge, inode1, inode2));
        } else {
          iedge = ~(itr->second);
        }
        CHECK((iedge >= 0 and mEdges[iedge].mNode1ID == inode1) or (mEdges[~iedge].mNode2ID == inode1));
        faceEdges.push_back(iedge);
      }
      CHECK(faceEdges.size() == n);
      const FaceHash fhash(faceNodes.begin(), faceNodes.end());
      int iface;
      const map<FaceHash, unsigned>::const_iterator itr = faceHash2ID.find(fhash);
      if (itr == faceHash2ID.end()) {
        iface = faceHash2ID.size();
        faceHash2ID[fhash] = iface;
        const FaceZoneHash& zones = faceZones[fhash];
        CHECK(zones.first == newZoneID);
        mFaces.push_back(Face(*this, iface, zones.first, zones.second, faceEdges));
      } else {
        iface = itr->second;
        CHECK((int)positiveID(mFaces[iface].mZone1ID) == newZoneID or
              (int)positiveID(mFaces[iface].mZone2ID) == newZoneID);
        if (mFaces[iface].mZone1ID == ~newZoneID or
            mFaces[iface].mZone2ID == ~newZoneID) iface = ~iface;
      }
      zoneFaces.push_back(iface);
    }
    CHECK(zoneFaces.size() == zoneFaceNodes.size());
    mZones.push_back(Zone(*this, newZoneID, zoneFaces));
  }

  // Post-conditions.
  ENSURE(mNodes.size() == mNodePositions.size());
  ENSURE(mZones.size() == numOldZones + newCells.size());
}

}

//------------------------------------------------------------------------------
// Instantiate the generic mesh non-inlined methods.
//------------------------------------------------------------------------------
#include "Mesh.cc"
template class Spheral::Mesh<Spheral::Dim<3> >;
