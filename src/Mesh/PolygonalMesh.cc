//---------------------------------Spheral++----------------------------------//
// PolygonalMesh -- 2-D mesh class.
//
// Created by JMO, Tue Nov 16 14:18:20 PST 2010
//----------------------------------------------------------------------------//
#ifdef USE_POLYTOPE
#include "polytope/polytope.hh"
#endif

#include "Mesh.hh"
#include "Utilities/DBC.hh"
#include "Distributed/Communicator.hh"
#include "Utilities/Timer.hh"

#include <limits>
#include <set>
#include <sstream>
using std::vector;
using std::map;
using std::set;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

namespace {

#ifdef USE_POLYTOPE
//------------------------------------------------------------------------------
// Internal worker method with common code for building from a 2D polytope
// tessellation.
//------------------------------------------------------------------------------
void buildFromPolytope(polytope::Tessellation<2, double>& tessellation,
                       const Dim<2>::Vector& xmin,
                       const Dim<2>::Vector& xmax,
                       const std::vector<Dim<2>::Vector>& generators,
                       Mesh<Dim<2> >& mesh,
                       std::vector<Dim<2>::Vector>& mNodePositions,
                       Mesh<Dim<2> >::NodeContainer& mNodes,
                       Mesh<Dim<2> >::EdgeContainer& mEdges,
                       Mesh<Dim<2> >::FaceContainer& mFaces,
                       Mesh<Dim<2> >::ZoneContainer& mZones,
                       std::vector<unsigned>& mNeighborDomains,
                       std::vector<std::vector<unsigned> > mSharedNodes,
                       std::vector<std::vector<unsigned> > mSharedFaces) {

  typedef Dim<2>::Vector Vector;
  typedef Mesh<Dim<2> > MeshType;
  typedef MeshType::Node Node;
  typedef MeshType::Edge Edge;
  typedef MeshType::Face Face;
  typedef MeshType::Zone Zone;
  const unsigned UNSETID = MeshType::UNSETID;

  // Read out and construct the Mesh node positions.
  TIME_FUNCTION
  int i, j, k, igen, jgen;
  const unsigned numGens = generators.size();
  const unsigned numNodes = tessellation.nodes.size()/2;
  mNodePositions.reserve(numNodes);
  for (i = 0; i != numNodes; ++i) {
    mNodePositions.push_back(Vector(max(xmin.x(), min(xmax.x(), tessellation.nodes[2*i])),
                                    max(xmin.y(), min(xmax.y(), tessellation.nodes[2*i+1]))));
  }

  // Build the edges and faces, which are degnerate in 2D.
  // Simultaneously build up the node->zone connectivity.
  const unsigned numEdges = tessellation.faces.size();
  unsigned inode, jnode;
  map<unsigned, set<unsigned> > nodeZones;
  for (i = 0; i != numEdges; ++i) {
    CHECK(tessellation.faces[i].size() == 2);
    CHECK(tessellation.faceCells[i].size() == 1 or
          tessellation.faceCells[i].size() == 2);
    inode = tessellation.faces[i][0];
    jnode = tessellation.faces[i][1];
    igen = tessellation.faceCells[i][0];
    jgen = (tessellation.faceCells[i].size() == 2 ? 
            tessellation.faceCells[i][1] :
            igen < 0 ?
            UNSETID:
            ~UNSETID);

    // Do we need to flip this face?
    if (igen < 0 and jgen == UNSETID) {
      igen = ~igen;
      jgen = ~jgen;
      swap(inode, jnode);
      vector<int>::iterator itr = find(tessellation.cells[igen].begin(), tessellation.cells[igen].end(), ~i);
      CHECK(itr != tessellation.cells[igen].end());
      *itr = i;
    }

    mEdges.push_back(Edge(mesh, i, inode, jnode));
    mFaces.push_back(Face(mesh, i, igen, jgen, vector<unsigned>(1, i)));

    for (const int j: tessellation.faceCells[i]) {
      nodeZones[inode].insert(MeshType::positiveID(j));
      nodeZones[jnode].insert(MeshType::positiveID(j));
    }
  }
  CHECK(mEdges.size() == numEdges);
  CHECK(mFaces.size() == numEdges);
  CHECK(nodeZones.size() == numNodes);

  // Construct the nodes.
  mNodes.reserve(numNodes);
  for (i = 0; i != numNodes; ++i) mNodes.push_back(Node(mesh, i, vector<unsigned>(nodeZones[i].begin(),
                                                                                  nodeZones[i].end())));
  CHECK(mNodes.size() == numNodes);

  // Construct the zones.
  mZones.reserve(numGens);
  for (i = 0; i != numGens; ++i) mZones.push_back(Zone(mesh, i, tessellation.cells[i]));
  CHECK(mZones.size() == numGens);

  // Copy the parallel info.
  mNeighborDomains = tessellation.neighborDomains;
  mSharedNodes = tessellation.sharedNodes;
  mSharedFaces = tessellation.sharedFaces;
  // cerr << "Assigned neighbor domain info : " << mNeighborDomains.size() << " " << mSharedNodes.size() << " " << mSharedFaces.size() << endl;

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    // Make sure any faces on the surface are pointing out of the mesh.
    for (const Face& face: mFaces) {
      ENSURE(face.zone1ID() != UNSETID);
      ENSURE(face.zone2ID() != UNSETID);
      ENSURE(MeshType::positiveID(face.zone2ID()) != UNSETID or face.zone1ID() >= 0);
      // ENSURE2(MeshType::positiveID(face.zone2ID()) != UNSETID or
      //         face.unitNormal().dot(mZones[face.zone1ID()].position() - face.position()) < 0.0,
      //         "Something amiss at the surface of the mesh : "
      //         << face.zone1ID() << " " << face.zone2ID() << " : " 
      //         << face.unitNormal() << " dot (" << mZones[face.zone1ID()].position() << " - " << face.position() << ") = "
      //         << face.unitNormal().dot(mZones[face.zone1ID()].position() - face.position()));
    }
    ENSURE2(mesh.valid() == "", mesh.valid());
    ENSURE2(mesh.validDomainInfo(xmin, xmax, false) == "", mesh.validDomainInfo(xmin, xmax, false));
  }
  END_CONTRACT_SCOPE

}
#endif

}

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<> const unsigned Mesh<Dim<2> >::minFacesPerZone = 3;
template<> const unsigned Mesh<Dim<2> >::minEdgesPerZone = 3;
template<> const unsigned Mesh<Dim<2> >::minNodesPerZone = 3;
template<> const unsigned Mesh<Dim<2> >::minEdgesPerFace = 1;
template<> const unsigned Mesh<Dim<2> >::minNodesPerFace = 2;

//------------------------------------------------------------------------------
// Mesh::reconstructInternal (xmin, xmax)
//------------------------------------------------------------------------------
template<>
void
Mesh<Dim<2> >::
reconstructInternal(const vector<Dim<2>::Vector>& generators,
                    const Dim<2>::Vector& xmin,
                    const Dim<2>::Vector& xmax) {
  CONTRACT_VAR(generators);
  CONTRACT_VAR(xmin);
  CONTRACT_VAR(xmax);

#ifdef USE_POLYTOPE
  // Some useful typedefs.
  typedef Dim<2> Dimension;

  // Is there anything to do?
  if (generators.size() == 0) return;

  // The tolerance on which we consider positions to be degenerate.
  const double xtol = 1.0e-12*(xmax - xmin).maxElement();
  const double xtol2 = xtol*xtol;

  // Pre-conditions.
  int i, j, k, igen, numGens = generators.size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(xmin.x() < xmax.x() and
            xmin.y() < xmax.y());
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

  // Copy the generator positions to a polytope style flat array.
  vector<double> gens;
  gens.reserve(2*generators.size());
  for (igen = 0; igen != numGens; ++igen) {
    gens.push_back(generators[igen].x());
    gens.push_back(generators[igen].y());
  }
  CHECK(gens.size() == 2*numGens);

  // Do the polytope tessellation.
  TIME_FUNCTION
  polytope::Tessellation<2, double> tessellation;
  {
#if 0     //#ifdef SPHERAL_ENABLE_MPI                                // FIXME when parallel polytope working again!
    polytope::DistributedTessellator<2, double> tessellator
#if defined USE_TRIANGLE && ( USE_TRIANGLE>0 )
      (new polytope::TriangleTessellator<double>(),
#else
      (new polytope::BoostTessellator<double>(),
#endif
       true,     // Manage memory for serial tessellator
       true);    // Build parallel connectivity
#else
#if defined USE_TRIANGLE && ( USE_TRIANGLE>0 )
    polytope::TriangleTessellator<double> tessellator;
#else
    polytope::BoostTessellator<double> tessellator;
#endif
    tessellator.tessellate(gens, const_cast<double*>(xmin.begin()), const_cast<double*>(xmax.begin()), tessellation);
#endif
    tessellator.tessellate(gens, const_cast<double*>(xmin.begin()), const_cast<double*>(xmax.begin()), tessellation);
  }
  CHECK(tessellation.cells.size() == numGens);

  // // Blago!
  // {
  //   vector<double> index(tessellation.cells.size());
  //   for (int i = 0; i < tessellation.cells.size(); ++i) index[i] = double(i);
  //   map<string, double*> fields;
  //   fields["cell_index"] = &index[0];
  //   polytope::SiloWriter<2, double>::write(tessellation, fields, "polygonal_blago");
  // }
  // // Blago!

  // Dispatch the work to our common 2D coding.
  buildFromPolytope(tessellation,
                    xmin,
                    xmax,
                    generators,
                    *this, 
                    mNodePositions,
                    mNodes,
                    mEdges,
                    mFaces,
                    mZones,
                    mNeighborDomains,
                    mSharedNodes,
                    mSharedFaces);
#endif
}

//------------------------------------------------------------------------------
// Mesh::reconstructInternal (FacetedVolume)
//------------------------------------------------------------------------------
template<>
void
Mesh<Dim<2> >::
reconstructInternal(const vector<Dim<2>::Vector>& generators,
                    const Dim<2>::FacetedVolume& boundary) {
  CONTRACT_VAR(generators);
  CONTRACT_VAR(boundary);
#ifdef USE_POLYTOPE

  // Some useful typedefs.
  typedef Dim<2> Dimension;

  // Is there anything to do?
  if (generators.size() == 0) return;

  // The tolerance on which we consider positions to be degenerate.
  const Vector xmin = boundary.xmin(), xmax = boundary.xmax();
  const double xtol = 1.0e-12*(xmax - xmin).maxElement();
  const double xtol2 = xtol*xtol;

  // Pre-conditions.
  int i, j, k, igen, numGens = generators.size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(xmin.x() < xmax.x() and
            xmin.y() < xmax.y());
    for (igen = 0; igen != numGens; ++igen)
      REQUIRE2(boundary.contains(generators[igen]), "Generator out of bounds " << generators[igen]);
    for (i = 0; i < numGens - 1; ++i) {
      for (j = i + 1; j < numGens; ++j) {
        REQUIRE2((generators[i] - generators[j]).magnitude2() > xtol2, 
                 "Degenerate generator positions:  " << i << " " << j << " " << generators[i] << " " << generators[j]);
      }
    }
  }
  END_CONTRACT_SCOPE

  // Copy the generator positions to a polytope style flat array.
  vector<double> gens;
  gens.reserve(2*generators.size());
  for (igen = 0; igen != numGens; ++igen) {
    gens.push_back(generators[igen].x());
    gens.push_back(generators[igen].y());
  }
  CHECK(gens.size() == 2*numGens);

  // Create a polytope PLC from our faceted boundary.
  polytope::ReducedPLC<2, double> plcBoundary;
  {
    const vector<Vector>& vertices = boundary.vertices();
    plcBoundary.points.reserve(2*vertices.size());
    for (vector<Vector>::const_iterator itr = vertices.begin();
         itr != vertices.end();
         ++itr) {
      plcBoundary.points.push_back(itr->x());
      plcBoundary.points.push_back(itr->y());
    }
    const vector<vector<unsigned> > facetVertices = boundary.facetVertices();
    plcBoundary.facets = vector<vector<int> >(facetVertices.size());
    for (i = 0; i != facetVertices.size(); ++i) {
      CHECK(facetVertices[i].size() == 2);
      for (j = 0; j != 2; ++j) {
        plcBoundary.facets[i].push_back(facetVertices[i][j]);
      }
    }
  }

  // Do the polytope tessellation.
  TIME_FUNCTION
  polytope::Tessellation<2, double> tessellation;
  {
#if 0   //  #ifdef SPHERAL_ENABLE_MPI                                    // FIXME when polytope Distributed fixed
    polytope::DistributedTessellator<2, double> tessellator
#if defined USE_TRIANGLE && ( USE_TRIANGLE>0 )
      (new polytope::TriangleTessellator<double>(),
#else
      (new polytope::BoostTessellator<double>(),
#endif
                                                            true,     // Manage memory for serial tessellator
                                                            true);    // Build parallel connectivity
#else
#if defined USE_TRIANGLE && ( USE_TRIANGLE>0 )
    polytope::TriangleTessellator<double> tessellator;
#else
    polytope::BoostTessellator<double> tessellator;
#endif
    tessellator.tessellate(gens, const_cast<double*>(xmin.begin()), const_cast<double*>(xmax.begin()), tessellation);
#endif
    tessellator.tessellate(gens, plcBoundary.points, plcBoundary, tessellation);
  }
  CHECK(tessellation.cells.size() == numGens);

  // // Blago!
  // {
  //   cerr << "Writing crap out..." << endl;
  //   vector<double> index(tessellation.cells.size());
  //   for (int i = 0; i < tessellation.cells.size(); ++i) index[i] = double(i);
  //   map<string, double*> fields, nullfields;
  //   fields["cell_index"] = &index[0];
  //   polytope::SiloWriter<2, double>::write(tessellation, 
  //                                          nullfields, nullfields, nullfields, 
  //                                          fields, "polygonal_blago");
  // }
  // // Blago!

  // Dispatch the work to our common 2D coding.
  buildFromPolytope(tessellation,
                    xmin,
                    xmax,
                    generators,
                    *this, 
                    mNodePositions,
                    mNodes,
                    mEdges,
                    mFaces,
                    mZones,
                    mNeighborDomains,
                    mSharedNodes,
                    mSharedFaces);

#endif
}

//------------------------------------------------------------------------------
// Compute the bounding surface of the mesh.
//------------------------------------------------------------------------------
template<>
Dim<2>::FacetedVolume
Mesh<Dim<2> >::
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
  int id2;
  unsigned i, j, iglobal, jglobal;
  map<unsigned, Vector> globalVertexPositions;
  vector<vector<unsigned> > facetIndices;
  for (const Face& face: mFaces) {
    CHECK(face.mNodeIDs.size() == 2);
    i = face.mNodeIDs[0];
    j = face.mNodeIDs[1];
    //id1 = face.zone1ID();
    id2 = face.zone2ID();
    if (positiveID(id2) == UNSETID and
        (sharedNodes.find(i) == sharedNodes.end() or
         sharedNodes.find(j) == sharedNodes.end())) {
      iglobal = local2globalIDs[i];
      jglobal = local2globalIDs[j];
      globalVertexPositions[iglobal] = mNodePositions[i];
      globalVertexPositions[jglobal] = mNodePositions[j];
      vector<unsigned> ids;
      ids.push_back(iglobal);
      ids.push_back(jglobal);
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
          CHECK2(otherIndices.size() == 2, "Bad size: " << otherIndices.size());
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
    CHECK(indices.size() == 2);
    CHECK(global2vertexID.find(indices[0]) != global2vertexID.end());
    CHECK(global2vertexID.find(indices[1]) != global2vertexID.end());
    indices[0] = global2vertexID[indices[0]];
    indices[1] = global2vertexID[indices[1]];
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    for (const vector<unsigned>& indices: facetIndices) {
      CONTRACT_VAR(indices);
      ENSURE(indices.size() == 2);
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
Mesh<Dim<2> >::
createNewMeshElements(const vector<vector<vector<unsigned> > >& newCells) {

  typedef pair<unsigned, unsigned> EdgeHash;
  typedef pair<int, int> FaceZoneHash;

  // Pre-conditions.
  REQUIRE(mNodes.size() <= mNodePositions.size());
  BEGIN_CONTRACT_SCOPE
  {
    for (const vector<vector<unsigned> >& cellFaces: newCells) {
      REQUIRE(cellFaces.size() >= minFacesPerZone);
      for (const vector<unsigned>& faceNodes: cellFaces) {
        REQUIRE(faceNodes.size() >= minNodesPerFace);
        for (unsigned inode: faceNodes) {
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
  map<EdgeHash, FaceZoneHash> faceZones;
  for (unsigned iface = 0; iface != numOldFaces; ++iface) {
    CHECK(mFaces[iface].mNodeIDs.size() == 2);
    faceZones[hashEdge(mFaces[iface].mNodeIDs[0], mFaces[iface].mNodeIDs[1])] = FaceZoneHash(mFaces[iface].mZone1ID, mFaces[iface].mZone2ID);
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
      CHECK(faceNodes.size() == 2);
      // cerr << "Adding face with nodes : ";
      // copy(faceNodes.begin(), faceNodes.end(), ostream_iterator<unsigned>(cerr, " "));
      const EdgeHash fhash = hashEdge(faceNodes[0], faceNodes[1]);
      map<EdgeHash, FaceZoneHash>::iterator faceItr = faceZones.find(fhash);
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
  for (map<EdgeHash, FaceZoneHash>::const_iterator faceItr = faceZones.begin();
       faceItr != faceZones.end();
       ++faceItr) {
    const EdgeHash& nodes = faceItr->first;
    const unsigned zone1 = positiveID(faceItr->second.first);
    const unsigned zone2 = positiveID(faceItr->second.second);
    nodeZones[nodes.first].insert(zone1);
    nodeZones[nodes.first].insert(zone2);
    nodeZones[nodes.second].insert(zone1);
    nodeZones[nodes.second].insert(zone2);
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
  // Since this is also the face hash->faceID mapping, 
  // we simultaneously update the face->zone connectivity.
  map<EdgeHash, unsigned> edgeHash2ID;
  for (unsigned iedge = 0; iedge != numOldEdges; ++iedge) {
    EdgeHash ehash = hashEdge(mEdges[iedge].mNode1ID, mEdges[iedge].mNode2ID);
    edgeHash2ID[ehash] = iedge;
    const FaceZoneHash& zones = faceZones[ehash];
    mFaces[iedge].mZone1ID = zones.first;
    mFaces[iedge].mZone2ID = zones.second;
  }

  // Create any new edges, faces, and zones.
  for (unsigned i = 0; i != newCells.size(); ++i) {
    const int newZoneID = numOldZones + i;
    const vector<vector<unsigned> >& zoneFaceNodes = newCells[i];
    vector<int> zoneFaces;
    for (const vector<unsigned>& faceNodes: zoneFaceNodes) {
      CHECK(faceNodes.size() == 2);
      const unsigned inode1 = faceNodes[0];
      const unsigned inode2 = faceNodes[1];
      const EdgeHash ehash = hashEdge(inode1, inode2);
      int iedge;
      const map<EdgeHash, unsigned>::const_iterator itr = edgeHash2ID.find(ehash);
      if (itr == edgeHash2ID.end()) {
        iedge = edgeHash2ID.size();
        edgeHash2ID[ehash] = iedge;
        mEdges.push_back(Edge(*this, iedge, inode1, inode2));
        const FaceZoneHash& zones = faceZones[ehash];
        CHECK(zones.first == newZoneID);
        mFaces.push_back(Face(*this, iedge, zones.first, zones.second, vector<unsigned>(1, iedge)));
      } else {
        iedge = ~(itr->second);
      }
      CHECK((iedge >= 0 and mEdges[iedge].mNode1ID == inode1) or (mEdges[~iedge].mNode2ID == inode1));
      zoneFaces.push_back(iedge);
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
template class Spheral::Mesh<Spheral::Dim<2> >;
