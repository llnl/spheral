//---------------------------------Spheral++----------------------------------//
// Mesh
//
// Created by JMO, Tue Oct 12 23:07:22 PDT 2010
//----------------------------------------------------------------------------//
#include <map>
#include <unordered_map>
#include <functional>
#include "boost/bimap.hpp"

using std::unordered_map;

#include "MeshConstructionUtilities.hh"
#include "Utilities/removeElements.hh"
#include "Utilities/DBC.hh"
#include "Distributed/allReduce.hh"
#include "Utilities/boundingBox.hh"
#include "Distributed/Communicator.hh"
#include "NodeList/NodeList.hh"

#ifdef SPHERAL_ENABLE_MPI
#include <mpi.h>
#include "Utilities/packElement.hh"
#endif

#include <limits>
#include <numeric>
#include <list>
using std::vector;
using std::set;
using std::list;
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

namespace { // anonymous

#ifdef SPHERAL_ENABLE_MPI
//------------------------------------------------------------------------------
// Look for anyone who has a non-empty string.
//------------------------------------------------------------------------------
string
reduceToMaxString(const string& x,
                  const unsigned rank,
                  const unsigned numDomains) {
  unsigned badRank = allReduce((x.size() == 0 ? numDomains : rank),
                               SPHERAL_OP_MIN);
  if (badRank == numDomains) {
    return "";
  } else {
    unsigned size = x.size();
    MPI_Bcast(&size, 1, MPI_UNSIGNED, badRank, Communicator::communicator());
    vector<char> result(x.begin(), x.end());
    result.resize(size);
    MPI_Bcast(&result.front(), size, MPI_CHAR, badRank, Communicator::communicator());
    return string(result.begin(), result.end());
  }
}
#endif

}           // anonymous

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
Mesh<Dimension>::
Mesh():
  mNodePositions(),
  mNodes(),
  mEdges(),
  mFaces(),
  mZones(),
  mNeighborDomains(),
  mSharedNodes(),
  mSharedFaces(),
  mNodeListNameOffsets(),
  mNodeListIndexOffsets() {
}

//------------------------------------------------------------------------------
// Mesh::Mesh(generators, xmin, xmax)
//------------------------------------------------------------------------------
template<typename Dimension>
Mesh<Dimension>::
Mesh(const vector<typename Dimension::Vector>& generators,
     const typename Dimension::Vector& xmin,
     const typename Dimension::Vector& xmax):
  mNodePositions(),
  mNodes(),
  mEdges(),
  mFaces(),
  mZones(),
  mNeighborDomains(),
  mSharedNodes(),
  mSharedFaces(),
  mNodeListNameOffsets(),
  mNodeListIndexOffsets() {
  this->reconstruct(generators, xmin, xmax);
}

//------------------------------------------------------------------------------
// Mesh::Mesh(generators, polygon)
//------------------------------------------------------------------------------
template<typename Dimension>
Mesh<Dimension>::
Mesh(const vector<typename Dimension::Vector>& generators,
     const typename Dimension::FacetedVolume& boundary):
  mNodePositions(),
  mNodes(),
  mEdges(),
  mFaces(),
  mZones(),
  mNeighborDomains(),
  mSharedNodes(),
  mSharedFaces(),
  mNodeListNameOffsets(),
  mNodeListIndexOffsets() {
  this->reconstruct(generators, boundary);
}

//------------------------------------------------------------------------------
// The explict constructor where the user directly passes in all the 
// information.
//------------------------------------------------------------------------------
template<typename Dimension>
Mesh<Dimension>::
Mesh(const vector<Vector>& nodePositions,
     const vector<vector<unsigned> >& edgeNodes,
     const vector<vector<unsigned> >& faceEdges,
     const vector<vector<int> >& zoneFaces):
  mNodePositions(nodePositions),
  mNodes(),
  mEdges(),
  mFaces(),
  mZones(),
  mNeighborDomains(),
  mSharedNodes(),
  mSharedFaces(),
  mNodeListNameOffsets(),
  mNodeListIndexOffsets() {

  // Invert the zoneFace structure.
  vector<vector<int> > faceZones(faceEdges.size());
  for (auto izone = 0u; izone != zoneFaces.size(); ++izone) {
    const vector<int>& zf = zoneFaces[izone];
    for (unsigned j = 0; j != zf.size(); ++j) {
      if (zf[j] < 0) {
        const int k = ~zf[j];
        VERIFY(k < (int)faceZones.size());
        faceZones[k].push_back(~izone);
      } else {
        VERIFY(zf[j] < (int)faceZones.size());
        faceZones[zf[j]].push_back(izone);
      }
    }
  }

  // Back out the set of zones for each node.
  vector<set<unsigned> > nodeZones(nodePositions.size());
  for (unsigned izone = 0; izone != zoneFaces.size(); ++izone) {
    const vector<int>& zf = zoneFaces[izone];
    for (unsigned i = 0; i != zf.size(); ++i) {
      const unsigned iface = positiveID(zf[i]);
      CHECK(iface < faceEdges.size());
      const vector<unsigned>& fe = faceEdges[iface];
      for (unsigned j = 0; j != fe.size(); ++j) {
        const unsigned iedge = fe[j];
        CHECK(iedge < edgeNodes.size());
        const vector<unsigned>& en = edgeNodes[iedge];
        CHECK(en.size() == 2);
        CHECK(en[0] < nodeZones.size());
        CHECK(en[1] < nodeZones.size());
        nodeZones[en[0]].insert(izone);
        nodeZones[en[1]].insert(izone);
      }
    }
  }

  // Construct the nodes.
  for (unsigned inode = 0; inode != nodePositions.size(); ++inode) {
    mNodes.push_back(Node(*this, inode, vector<unsigned>(nodeZones[inode].begin(), nodeZones[inode].end())));
  }

  // Construct the edges.
  for (unsigned iedge = 0; iedge != edgeNodes.size(); ++iedge) {
    VERIFY(edgeNodes[iedge].size() == 2);
    mEdges.push_back(Edge(*this, iedge, edgeNodes[iedge][0], edgeNodes[iedge][1]));
  }

  // Construct the faces.
  for (unsigned iface = 0; iface != faceEdges.size(); ++iface) {
    VERIFY(faceZones[iface].size() == 1 or faceZones[iface].size() == 2);
    if (faceZones[iface].size() == 1) {
      mFaces.push_back(Face(*this, iface, faceZones[iface][0], UNSETID, faceEdges[iface]));
    } else {
      mFaces.push_back(Face(*this, iface, faceZones[iface][0], faceZones[iface][1], faceEdges[iface]));
    }
  }

  // Construct the zones.
  for (unsigned izone = 0; izone != zoneFaces.size(); ++izone) {
    mZones.push_back(Zone(*this, izone, zoneFaces[izone]));
  }

  // Post-conditions.
  ENSURE(mNodePositions.size() == nodePositions.size());
  ENSURE(mNodes.size() == nodePositions.size());
  ENSURE(mEdges.size() == edgeNodes.size());
  ENSURE(mFaces.size() == faceEdges.size());
  ENSURE(mZones.size() == zoneFaces.size());
}

//------------------------------------------------------------------------------
// Copy
//------------------------------------------------------------------------------
template<typename Dimension>
Mesh<Dimension>::
Mesh(const Mesh<Dimension>& rhs) {
  *this = rhs;
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<typename Dimension>
Mesh<Dimension>&
Mesh<Dimension>::
operator=(const Mesh<Dimension>& rhs) {
  if (this != &rhs) {
    mNodePositions = rhs.mNodePositions;
    mNodes = rhs.mNodes;
    mEdges = rhs.mEdges;
    mFaces = rhs.mFaces;
    mZones = rhs.mZones;
    mNeighborDomains = rhs.mNeighborDomains;
    mSharedNodes = rhs.mSharedNodes;
    mSharedFaces = rhs.mSharedFaces;
    mNodeListNameOffsets = rhs.mNodeListNameOffsets;
    mNodeListIndexOffsets = rhs.mNodeListIndexOffsets;

    // Set the mesh pointers appropriately.
    for (typename vector<Node>::iterator itr = mNodes.begin();
         itr != mNodes.end();
         ++itr) itr->mMeshPtr = this;
    for (typename vector<Edge>::iterator itr = mEdges.begin();
         itr != mEdges.end();
         ++itr) itr->mMeshPtr = this;
    for (typename vector<Face>::iterator itr = mFaces.begin();
         itr != mFaces.end();
         ++itr) itr->mMeshPtr = this;
    for (typename vector<Zone>::iterator itr = mZones.begin();
         itr != mZones.end();
         ++itr) itr->mMeshPtr = this;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Mesh::clear
//------------------------------------------------------------------------------
template<typename Dimension>
void
Mesh<Dimension>::
clear() {
  mNodePositions =        vector<Vector>();
  mNodes =                NodeContainer();
  mEdges =                EdgeContainer();
  mFaces =                FaceContainer();
  mZones =                ZoneContainer();
  mNeighborDomains =      vector<unsigned>();
  mSharedNodes =          vector<vector<unsigned> >();
  mSharedFaces =          vector<vector<unsigned> >();
  mNodeListNameOffsets =  map<string, unsigned>();
  mNodeListIndexOffsets = vector<unsigned>();
}

//------------------------------------------------------------------------------
// Mesh::reconstruct(generators, xmin, xmax)
//------------------------------------------------------------------------------
template<typename Dimension>
void
Mesh<Dimension>::
reconstruct(const vector<typename Dimension::Vector>& generators,
            const typename Dimension::Vector& xmin,
            const typename Dimension::Vector& xmax) {
  this->clear();
  this->reconstructInternal(generators, xmin, xmax);
}

//------------------------------------------------------------------------------
// Mesh::reconstruct(generators, polygon)
//------------------------------------------------------------------------------
template<typename Dimension>
void
Mesh<Dimension>::
reconstruct(const vector<typename Dimension::Vector>& generators,
            const typename Dimension::FacetedVolume& boundary) {
  this->clear();
  this->reconstructInternal(generators, boundary);
}

//------------------------------------------------------------------------------
// Mesh::removeZonesByMask
//------------------------------------------------------------------------------
template<typename Dimension>
void
Mesh<Dimension>::
removeZonesByMask(const vector<unsigned>& zoneMask) {

  // Pre-conditions.
  VERIFY2(zoneMask.size() == mZones.size(), "Mask wrong size:  " << zoneMask.size() << " != " << mZones.size());

  // Construct masks indicating which nodes, edges, and faces we need to keep
  // or remove.
  vector<unsigned> nodeMask(mNodes.size(), 0);
  vector<unsigned> edgeMask(mEdges.size(), 0);
  vector<unsigned> faceMask(mFaces.size(), 0);
  for (size_t izone = 0; izone != zoneMask.size(); ++izone) {
    if (zoneMask[izone] == 1) {
      const Zone& zone = mZones[izone];
      const vector<unsigned>& nodeIDs = zone.nodeIDs();
      const vector<unsigned>& edgeIDs = zone.edgeIDs();
      const vector<int>& faceIDs = zone.faceIDs();
      for (vector<unsigned>::const_iterator itr = nodeIDs.begin();
           itr != nodeIDs.end();
           ++itr) {
        CHECK(*itr < mNodes.size());
        nodeMask[*itr] = 1;
      }
      for (vector<unsigned>::const_iterator itr = edgeIDs.begin();
           itr != edgeIDs.end();
           ++itr) {
        CHECK(*itr < mEdges.size());
        edgeMask[*itr] = 1;
      }
      for (vector<int>::const_iterator itr = faceIDs.begin();
           itr != faceIDs.end();
           ++itr) {
        const int fid = positiveID(*itr);
        CHECK(fid >= 0 and fid < (int)mFaces.size());
        faceMask[fid] = 1;
      }
    }
  }
  
  // We need the new IDs for all elements.
  const vector<unsigned> newNodeIDs = this->recomputeIDs(nodeMask);
  const vector<unsigned> newEdgeIDs = this->recomputeIDs(edgeMask);
  const vector<unsigned> newFaceIDs = this->recomputeIDs(faceMask);
  const vector<unsigned> newZoneIDs = this->recomputeIDs(zoneMask);
  CHECK2(newNodeIDs.size() == mNodes.size(), newNodeIDs.size() << " " << mNodes.size());
  CHECK2(newEdgeIDs.size() == mEdges.size(), newEdgeIDs.size() << " " << mEdges.size());
  CHECK2(newFaceIDs.size() == mFaces.size(), newFaceIDs.size() << " " << mFaces.size());
  CHECK2(newZoneIDs.size() == mZones.size(), newZoneIDs.size() << " " << mZones.size());

  // Update the IDs of nodes and their internal data.
  for (unsigned i = 0; i != mNodes.size(); ++i) {
    Node& node = mNodes[i];
    node.mID = newNodeIDs[i];
    this->reassignIDs(node.mZoneIDs, newZoneIDs);
  }
  
  // Update the IDs of edges and their internal data.
  for (unsigned i = 0; i != mEdges.size(); ++i) {
    Edge& edge = mEdges[i];
    edge.mID = newEdgeIDs[i];
    edge.mNode1ID = newNodeIDs[edge.mNode1ID];
    edge.mNode2ID = newNodeIDs[edge.mNode2ID];
  }

  // Update the IDs of faces and their internal data.
  int z1ID, z2ID, s1, s2;
  for (unsigned i = 0; i != mFaces.size(); ++i) {
    Face& face = mFaces[i];
    face.mID = newFaceIDs[i];
    s1 = isgn(face.mZone1ID);
    s2 = isgn(face.mZone2ID);
    z1ID = positiveID(face.mZone1ID);
    z2ID = positiveID(face.mZone2ID);
    if (z1ID != UNSETID) face.mZone1ID = (s1 == 1 ? newZoneIDs[z1ID] : ~newZoneIDs[z1ID]);
    if (z2ID != UNSETID) face.mZone2ID = (s2 == 1 ? newZoneIDs[z2ID] : ~newZoneIDs[z2ID]);
    this->reassignIDs(face.mNodeIDs, newNodeIDs);
    this->reassignIDs(face.mEdgeIDs, newEdgeIDs);
  }

  // Update the IDs of the zones and their internal data.
  for (unsigned i = 0; i != mZones.size(); ++i) {
    Zone& zone = mZones[i];
    zone.mID = newZoneIDs[i];
    this->reassignIDs(zone.mNodeIDs, newNodeIDs);
    this->reassignIDs(zone.mEdgeIDs, newEdgeIDs);
    this->reassignIDs(zone.mFaceIDs, newFaceIDs);
  }

  // Erase the nodes.
  {
    vector<unsigned> kill;
    for (unsigned i = 0; i != nodeMask.size(); ++i) {
      if (nodeMask[i] == 0) kill.push_back(i);
    }
    removeElements(mNodes, kill);
    removeElements(mNodePositions, kill);
  }

  // Erase the edges.
  {
    vector<unsigned> kill;
    for (unsigned i = 0; i != edgeMask.size(); ++i) {
      if (edgeMask[i] == 0) kill.push_back(i);
    }
    removeElements(mEdges, kill);
  }

  // Erase the faces.
  {
    vector<unsigned> kill;
    for (unsigned i = 0; i != faceMask.size(); ++i) {
      if (faceMask[i] == 0) kill.push_back(i);
    }
    removeElements(mFaces, kill);
  }

  // Erase the zones.
  {
    vector<unsigned> kill;
    for (unsigned i = 0; i != zoneMask.size(); ++i) {
      if (zoneMask[i] == 0) kill.push_back(i);
    }
    removeElements(mZones, kill);
  }

  // Patch up the parallel information.
  vector<unsigned> killDomains;
  for (unsigned idomain = 0; idomain != mNeighborDomains.size(); ++idomain) {

    // Fix up the shared nodes.
    {
      vector<unsigned> kill;
      for (unsigned i = 0; i != mSharedNodes[idomain].size(); ++i) {
        if (nodeMask[mSharedNodes[idomain][i]] == 0) kill.push_back(i);
      }
      removeElements(mSharedNodes[idomain], kill);
      this->reassignIDs(mSharedNodes[idomain], newNodeIDs);
    }

    // Fix up the shared faces.
    {
      vector<unsigned> kill;
      for (unsigned i = 0; i != mSharedFaces[idomain].size(); ++i) {
        if (faceMask[mSharedFaces[idomain][i]] == 0) kill.push_back(i);
      }
      removeElements(mSharedFaces[idomain], kill);
      this->reassignIDs(mSharedFaces[idomain], newFaceIDs);
    }

    // Is there anything left or this domain?
    CHECK(mSharedFaces[idomain].size() <= mSharedNodes[idomain].size());
    if (mSharedNodes.size() == 0) killDomains.push_back(idomain);
  }
  removeElements(mNeighborDomains, killDomains);
  removeElements(mSharedNodes, killDomains);
  removeElements(mSharedFaces, killDomains);

  // // Any pre-existing parallel info is now invalid.
  // if (allReduce(mNeighborDomains.size(), SPHERAL_OP_MAX) > 0) {
  //   mNeighborDomains = vector<unsigned>();
  //   mSharedNodes = vector<vector<unsigned> >();
  //   mSharedFaces = vector<vector<unsigned> >();
  //   this->generateDomainInfo();
  // }

  // That's it.
  BEGIN_CONTRACT_SCOPE
  {
    Vector xmin, xmax;
    this->boundingBox(xmin, xmax);
    ENSURE2(this->validDomainInfo(xmin, xmax, false) == "", this->validDomainInfo(xmin, xmax, false));
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Mesh::cleanEdges
//------------------------------------------------------------------------------
template<typename Dimension>
void
Mesh<Dimension>::
cleanEdges(const double edgeTol) {

  // Pre-conditions.
  VERIFY2(edgeTol > 0.0, "Specify a positive (non-zero) edge tolerance!");

  // Iterate until all edges are above the threshold.
  bool meshClean = false;
  while (not meshClean) {
    meshClean = true;

    // Find the maximum edge length in each zone, and in turn find the maximum
    // zonal edge length for each edge of the zone.
    vector<double> maxZoneEdgeLength(mEdges.size(), 0.0);
    for (unsigned izone = 0; izone != mZones.size(); ++izone) {
      double maxZoneEdge = 0.0;
      vector<unsigned> edgeIDs = mZones[izone].edgeIDs();
      for (vector<unsigned>::const_iterator itr = edgeIDs.begin();
           itr != edgeIDs.end();
           ++itr) maxZoneEdge = max(maxZoneEdge, mEdges[*itr].length());
      for (vector<unsigned>::const_iterator itr = edgeIDs.begin();
           itr != edgeIDs.end();
           ++itr) maxZoneEdgeLength[*itr] = max(maxZoneEdgeLength[*itr], maxZoneEdge);
    }

    // Flag the edges we want to remove.  This implies we are also coalescing 
    // nodes.
    vector<unsigned> edgeMask(mEdges.size(), 1);
    vector<unsigned> nodeMask(mNodes.size(), 1);
    vector<unsigned> nodeMap(mNodes.size());
    for (unsigned i = 0; i != mNodes.size(); ++i) nodeMap[i] = i;
    unsigned n1, n2;
    for (unsigned iedge = 0; iedge != mEdges.size(); ++iedge) {
      n1 = mEdges[iedge].node1().ID();
      n2 = mEdges[iedge].node2().ID();
      if (mEdges[iedge].length() < edgeTol*maxZoneEdgeLength[iedge] and 
          nodeMask[n1] == 1 and nodeMask[n2] == 1) {
        meshClean = false;
        edgeMask[iedge] = 0;
        nodeMask[n1] = 2;
        nodeMask[n2] = 0;
        nodeMap[n2] = n1;
      }
    }
    replace_if(nodeMask.begin(), nodeMask.end(),
               [] (unsigned val) { return val == 2; },
               1);
    
    // Reassign the nodes we're renumbering before they get deleted.
    {
      // Edges.
      for (unsigned i = 0; i != mEdges.size(); ++i) {
        Edge& edge = mEdges[i];
        edge.mNode1ID = nodeMap[edge.mNode1ID];
        edge.mNode2ID = nodeMap[edge.mNode2ID];
      }

      // Faces.
      for (unsigned i = 0; i != mFaces.size(); ++i) {
        Face& face = mFaces[i];
        vector<unsigned> kill;
        for (unsigned j = 0; j != face.mNodeIDs.size(); ++j) {
          face.mNodeIDs[j] = nodeMap[face.mNodeIDs[j]];
        }
      }

      // // Zones.
      // for (unsigned i = 0; i != mZones.size(); ++i) {
      //   Zone& zone = mZones[i];
      //   for (unsigned j = 0; j != zone.mNodeIDs.size(); ++j) {
      //     zone.mNodeIDs[j] = nodeMap[zone.mNodeIDs[j]];
      //   }
      // }
    }

    // Check if there are any faces that need to be removed.
    vector<unsigned> faceMask(mFaces.size(), 1);
    unsigned numActiveEdges;
    for (unsigned iface = 0; iface != mFaces.size(); ++iface) {
      numActiveEdges = 0;
      const vector<unsigned>& edgeIDs = mFaces[iface].edgeIDs();
      for (vector<unsigned>::const_iterator itr = edgeIDs.begin();
           itr != edgeIDs.end();
           ++itr) numActiveEdges += edgeMask[*itr];
      faceMask[iface] = (numActiveEdges >= minEdgesPerFace ? 1 : 0);
    }

    // cout << "nodeMask: ";
    // copy(nodeMask.begin(), nodeMask.end(), ostream_iterator<unsigned>(cout, " "));
    // cout << endl
    //      << "edgeMask: ";
    // copy(edgeMask.begin(), edgeMask.end(), ostream_iterator<unsigned>(cout, " "));
    // cout << endl
    //      << "faceMask: ";
    // copy(faceMask.begin(), faceMask.end(), ostream_iterator<unsigned>(cout, " "));
    // cout << endl
    //      << "nodeMap: ";
    // copy(nodeMap.begin(), nodeMap.end(), ostream_iterator<unsigned>(cout, " "));
    // cout << endl;

    // Figure out the new numberings.
    const vector<unsigned> newNodeIDs = this->recomputeIDs(nodeMask);
    const vector<unsigned> newEdgeIDs = this->recomputeIDs(edgeMask);
    const vector<unsigned> newFaceIDs = this->recomputeIDs(faceMask);

    // cout << "newNodeIDs: ";
    // copy(newNodeIDs.begin(), newNodeIDs.end(), ostream_iterator<unsigned>(cout, " "));
    // cout << endl;
    // cout << "newEdgeIDs: ";
    // copy(newEdgeIDs.begin(), newEdgeIDs.end(), ostream_iterator<unsigned>(cout, " "));
    // cout << endl;
    // cout << "newFaceIDs: ";
    // copy(newFaceIDs.begin(), newFaceIDs.end(), ostream_iterator<unsigned>(cout, " "));
    // cout << endl;

    // Update the IDs of nodes and their internal data.
    {
      // cout << "Nodes:  " << endl;
      vector<unsigned> kill;
      for (unsigned i = 0; i != mNodes.size(); ++i) {
        if (newNodeIDs[i] == UNSETID) {
          kill.push_back(i);
        } else {
          Node& node = mNodes[i];
          node.mID = newNodeIDs[i];
          // cout << "    " << node.mID << " " << node.position() << endl;
        }
      }
      removeElements(mNodes, kill);
      removeElements(mNodePositions, kill);
    }
  
    // Update the IDs of edges and their internal data.
    {
      // cout << "Edges:  " << endl;
      vector<unsigned> kill;
      for (unsigned i = 0; i != mEdges.size(); ++i) {
        if (newEdgeIDs[i] == UNSETID) {
          kill.push_back(i);
        } else {
          Edge& edge = mEdges[i];
          edge.mID = newEdgeIDs[i];
          edge.mNode1ID = newNodeIDs[edge.mNode1ID];
          edge.mNode2ID = newNodeIDs[edge.mNode2ID];
          CHECK(edge.mNode1ID != UNSETID and
                edge.mNode2ID != UNSETID);
          // cout << "    " << edge.mID << " : "
          //      << edge.mNode1ID << " "
          //      << edge.mNode2ID << " "
          //      << endl;
        }
      }
      removeElements(mEdges, kill);
    }

    // Update the IDs of faces and their internal data.
    {
      // cout << "Faces:" << endl;
      vector<unsigned> kill;
      for (unsigned i = 0; i != mFaces.size(); ++i) {
        if (newFaceIDs[i] == UNSETID) {
          kill.push_back(i);
        } else {
          Face& face = mFaces[i];
          face.mID = newFaceIDs[i];
          this->reassignIDs(face.mNodeIDs, newNodeIDs);
          this->reassignIDs(face.mEdgeIDs, newEdgeIDs);
          this->removeUNSETIDs(face.mNodeIDs);
          this->removeUNSETIDs(face.mEdgeIDs);
          // cout << face.mID << " : ";
          // copy(face.mNodeIDs.begin(), face.mNodeIDs.end(), ostream_iterator<unsigned>(cout, " "));
          // cout << " : ";
          // copy(face.mEdgeIDs.begin(), face.mEdgeIDs.end(), ostream_iterator<unsigned>(cout, " "));
          // cout << endl;
        }
      }
      removeElements(mFaces, kill);
    }

    // Update the IDs of the zones and their internal data.
    {
      // cout << "Zones:  " << endl;
      vector<unsigned> kill;
      for (unsigned i = 0; i != mZones.size(); ++i) {
        Zone& zone = mZones[i];
        vector<int> faceIDs = zone.mFaceIDs;
        this->reassignIDs(faceIDs, newFaceIDs);
        this->removeUNSETIDs(faceIDs);
        zone = Zone(*this, i, faceIDs);
        // this->reassignIDs(zone.mEdgeIDs, newEdgeIDs);
        // this->reassignIDs(zone.mFaceIDs, newFaceIDs);
        // this->removeUNSETIDs(zone.mEdgeIDs);
        // this->removeUNSETIDs(zone.mFaceIDs);
        // zone.constructNodeIDs();
        // cout << zone.mID << " : ";
        // copy(zone.mNodeIDs.begin(), zone.mNodeIDs.end(), ostream_iterator<unsigned>(cout, " "));
        // cout << " : ";
        // copy(zone.mEdgeIDs.begin(), zone.mEdgeIDs.end(), ostream_iterator<unsigned>(cout, " "));
        // cout << " : ";
        // copy(zone.mFaceIDs.begin(), zone.mFaceIDs.end(), ostream_iterator<unsigned>(cout, " "));
        // cout << endl;
      }
    }
  }

  // Any pre-existing parallel info is now invalid, so just clear out the old
  // data.
  mNeighborDomains = vector<unsigned>();
  mSharedNodes = vector<vector<unsigned> >();
  mSharedFaces = vector<vector<unsigned> >();
}

//------------------------------------------------------------------------------
// Look up the NodeList offset and nodeID for the given zone.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Mesh<Dimension>::
lookupNodeListID(const unsigned zoneID, unsigned& nodeListi, unsigned& i) const {
  REQUIRE(zoneID < mZones.size());
  const std::vector<unsigned>::const_iterator itr = std::lower_bound(mNodeListIndexOffsets.begin(),
                                                                     mNodeListIndexOffsets.end(),
                                                                     zoneID);
  CHECK(itr != mNodeListIndexOffsets.end() and *itr >= zoneID);
  nodeListi = std::distance(mNodeListIndexOffsets.begin(), itr) - (*itr == zoneID ? 0 : 1);
  CHECK(zoneID >= mNodeListIndexOffsets[nodeListi]);
  i = zoneID - mNodeListIndexOffsets[nodeListi];
}

//------------------------------------------------------------------------------
// Mesh::generateDomainInfo
// Create the parallel domain info: i.e., which nodes on each domain are linked.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Mesh<Dimension>::
generateDomainInfo() {
  REQUIRE(mNodePositions.size() == numNodes());

  // This method is empty and a no-op unless we're building a parallel code!
#ifdef SPHERAL_ENABLE_MPI

  // Start out by determining the global extent of the mesh.
  Vector xmin, xmax;
  this->boundingBox(xmin, xmax);

  // Define the hashing scale.
  const double dxhash = (xmax - xmin).maxElement() / double(std::numeric_limits<KeyElement>::max());

  // Puff out the bounds a bit.  We do the all reduce just to ensure
  // bit perfect consistency across processors.
  Vector boxInv;
  for (unsigned i = 0; i != Dimension::nDim; ++i) {
    xmin(i) = allReduce(xmin(i) - dxhash, SPHERAL_OP_MIN);
    xmax(i) = allReduce(xmax(i) + dxhash, SPHERAL_OP_MAX);
    boxInv(i) = safeInv(xmax(i) - xmin(i));
  }

  // Hash the node positions.  We want these sorted by key as well
  // to make testing if a key is present fast.
  vector<Key> nodeHashes;
  std::map<Key, unsigned> key2nodeID;
  nodeHashes.reserve(numNodes());
  for (unsigned i = 0; i != numNodes(); ++i) {
    nodeHashes.push_back(hashPosition(mNodePositions[i], xmin, xmax, boxInv));
    key2nodeID[nodeHashes.back()] = i;
  }
  sort(nodeHashes.begin(), nodeHashes.end());
  CHECK2(nodeHashes.size() == numNodes(), "Bad sizes:  " << nodeHashes.size() << " " << numNodes());
  CHECK2(key2nodeID.size() == numNodes(), "Bad sizes:  " << key2nodeID.size() << " " << numNodes());

  // Puff out our domain positions a bit from their centroid to try and
  // ensure our intersection tests don't miss something.
  vector<Vector> hullPoints;
  {
    Vector centroid;
    for (unsigned i = 0; i != mNodePositions.size(); ++i) centroid += mNodePositions[i];
    CHECK(mNodePositions.size() > 0);
    centroid /= mNodePositions.size();
    hullPoints.reserve(mNodePositions.size());
    for (unsigned i = 0; i != mNodePositions.size(); ++i) 
      hullPoints.push_back(1.05*(mNodePositions[i] - centroid) + centroid);
  }

  // Generate convex hulls enclosing each domain.
  const unsigned numDomains = Process::getTotalNumberOfProcesses();
  const unsigned rank = Process::getRank();
  CHECK(rank < numDomains);
  vector<ConvexHull> domainHulls(numDomains);
  domainHulls[rank] = ConvexHull(hullPoints);

  // Globally exchange the convex hulls.  This might be a bottle neck!
  {
    vector<char> localBuffer;
    packElement(domainHulls[rank], localBuffer);
    for (auto sendProc = 0u; sendProc != numDomains; ++sendProc) {
      unsigned bufSize = localBuffer.size();
      MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, Communicator::communicator());
      vector<char> buffer = localBuffer;
      buffer.resize(bufSize);
      MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, Communicator::communicator());
      vector<char>::const_iterator itr = buffer.begin();
      unpackElement(domainHulls[sendProc], itr, buffer.end());
      CHECK(itr == buffer.end());
    }
  }

  // Check which bounding volume hulls intersect our own.  This represents the
  // set we'll potentially have to communicate with.
  vector<unsigned> potentialNeighborDomains;
  for (unsigned idomain = 0; idomain != numDomains; ++idomain) {
    if (idomain != rank and domainHulls[idomain].intersect(domainHulls[rank])) potentialNeighborDomains.push_back(idomain);
  }

  // Ensure consistency in the potential exchange pattern!
  BEGIN_CONTRACT_SCOPE
  {
    if (numDomains > 1) {
      for (int sendProc = 0; sendProc != (int)numDomains; ++sendProc) {
        unsigned num = potentialNeighborDomains.size();
        MPI_Bcast(&num, 1, MPI_UNSIGNED, sendProc, Communicator::communicator());
        vector<unsigned> otherNeighbors = potentialNeighborDomains;
        otherNeighbors.resize(num);
        MPI_Bcast(&otherNeighbors.front(), num, MPI_UNSIGNED, sendProc, Communicator::communicator());
        CHECK((binary_search(potentialNeighborDomains.begin(), potentialNeighborDomains.end(), sendProc) == true and
               binary_search(otherNeighbors.begin(), otherNeighbors.end(), rank) == true) or
              (binary_search(potentialNeighborDomains.begin(), potentialNeighborDomains.end(), sendProc) == false and
               binary_search(otherNeighbors.begin(), otherNeighbors.end(), rank) == false));
      }
    }
  }
  END_CONTRACT_SCOPE

  // Exchange the keys.
  vector<vector<Key> > neighborHashes;
  exchangeTuples(nodeHashes, potentialNeighborDomains, neighborHashes);

  // Determine which of our nodes are shared with each neighbor domain.
  const unsigned numNeighborDomains = potentialNeighborDomains.size();
  mNeighborDomains.reserve(numNeighborDomains);
  mSharedNodes.reserve(numNeighborDomains);
  for (unsigned k = 0; k != numNeighborDomains; ++k) {
    const unsigned otherProc = potentialNeighborDomains[k];
    const vector<Key>& otherHashes = neighborHashes[k];
    vector<unsigned> sharedIDs;
    sharedIDs.reserve(otherHashes.size());
    for (vector<Key>::const_iterator keyItr = nodeHashes.begin();
         keyItr != nodeHashes.end();
         ++keyItr) {
      if (binary_search(otherHashes.begin(), otherHashes.end(), *keyItr) == true) sharedIDs.push_back(key2nodeID[*keyItr]);
    }
    if (sharedIDs.size() > 0) {
      mNeighborDomains.push_back(otherProc);
      mSharedNodes.push_back(sharedIDs);
    }
  }
  CHECK(mNeighborDomains.size() == mSharedNodes.size());

  // Determine which process owns each shared node.  We use the convention that the process
  // with the lowest rank wins.
  std::unordered_map<unsigned, unsigned> nodeOwner;
  for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
    for (unsigned j = 0; j != mSharedNodes[k].size(); ++j) {
      CHECK(k < mSharedNodes.size() and j < mSharedNodes[k].size());
      nodeOwner[mSharedNodes[k][j]] = numDomains;
    }
  }
  for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
    const unsigned potentialOwner = min(rank, mNeighborDomains[k]);
    for (unsigned j = 0; j != mSharedNodes[k].size(); ++j) {
      CHECK(k < mSharedNodes.size() and j < mSharedNodes[k].size());
      const unsigned i = mSharedNodes[k][j];
      CHECK(nodeOwner.find(i) != nodeOwner.end());
      nodeOwner[i] = min(nodeOwner[i], potentialOwner);
    }
  }

  // Build the shared face info based on the nodes.
  mSharedFaces.resize(mNeighborDomains.size());
  for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
    set<unsigned> neighborNodes(mSharedNodes[k].begin(), mSharedNodes[k].end());
    for (typename vector<Face>::const_iterator faceItr = mFaces.begin();
         faceItr != mFaces.end();
         ++faceItr) {
      const vector<unsigned>& nodeIDs = faceItr->nodeIDs();
      vector<unsigned>::const_iterator itr = nodeIDs.begin();
      while (itr != nodeIDs.end() and 
             neighborNodes.find(*itr) != neighborNodes.end()) ++itr;
      if (itr == nodeIDs.end()) mSharedFaces[k].push_back(faceItr->ID());
    }
  }
  CHECK(mSharedFaces.size() == mNeighborDomains.size());

  // Check that the communicated info is consistent.
  CHECK2(validDomainInfo(xmin, xmax, false) == "", validDomainInfo(xmin, xmax, false));
        
  // Turns out we actually want to keep the duplicates nodes!
//   // Cull the lists of shared nodes so that we only talk with the owner about any individual
//   // node.
//   {
//     vector<unsigned> procsToKill;
//     for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
//       const unsigned otherProc = mNeighborDomains[k];
//       vector<unsigned> nodesToKill;
//       CHECK(k < mSharedNodes.size());
//       for (unsigned j = 0; j != mSharedNodes[k].size(); ++j) {
//         CHECK(k < mSharedNodes.size() and j < mSharedNodes[k].size());
//         const unsigned nodeID = mSharedNodes[k][j];
//         CHECK(nodeOwner.find(nodeID) != nodeOwner.end());
//         const unsigned owner = nodeOwner[nodeID];
//         CHECK(owner < numDomains);
//         if (owner != rank and owner != otherProc) nodesToKill.push_back(j);
//       }
//       if (nodesToKill.size() == mSharedNodes[k].size()) {
//         procsToKill.push_back(k);
//       } else {
//         removeElements(mSharedNodes[k], nodesToKill);
//       }
//     }
//     removeElements(mNeighborDomains, procsToKill);
//     removeElements(mSharedNodes, procsToKill);
//   }
//   CHECK(mNeighborDomains.size() == mSharedNodes.size());

//   // Check that the communicated info is consistent.
//   CHECK2(validDomainInfo(xmin, xmax, true) == "", validDomainInfo(xmin, xmax, true));

  // For the sake of bit-perfectness, exchange node positions such that every shared
  // node has the position it takes on the owner domain.
  {
    // Start by sending our info for all nodes we own.
    list<vector<char> > sendBuffers;
    list<unsigned> sendBufferSizes;
    vector<MPI_Request> sendRequests;
    sendRequests.reserve(2*mNeighborDomains.size());
    for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
      CHECK(k < mSharedNodes.size());
      const unsigned otherProc = mNeighborDomains[k];
      if (otherProc > rank) {
        vector<char> buf;
        for (unsigned j = 0; j != mSharedNodes[k].size(); ++j) {
          const unsigned i = mSharedNodes[k][j];
          CHECK(i < mNodePositions.size());
          CHECK(nodeOwner.find(i) != nodeOwner.end());
          if (nodeOwner[i] == rank) packElement(mNodePositions[i], buf);
        }
        sendBufferSizes.push_back(buf.size());
        sendBuffers.push_back(buf);
        sendRequests.push_back(MPI_Request());
        MPI_Isend(&sendBufferSizes.back(), 1, MPI_UNSIGNED, otherProc, 10, Communicator::communicator(), &sendRequests.back());
        if (buf.size() > 0) {
          sendRequests.push_back(MPI_Request());
          MPI_Isend(&(sendBuffers.back().front()), buf.size(), MPI_CHAR, otherProc, 11, Communicator::communicator(), &sendRequests.back());
        }
      }
    }
    CHECK(sendRequests.size() <= 2*mNeighborDomains.size());
    CHECK(sendBufferSizes.size() == sendBuffers.size());

    // Now get the postions from all domains sending to us.
    for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
      const unsigned otherProc = mNeighborDomains[k];
      if (otherProc < rank) {
        unsigned bufSize;
        MPI_Status stat1, stat2;
        MPI_Recv(&bufSize, 1, MPI_UNSIGNED, otherProc, 10, Communicator::communicator(), &stat1);
        if (bufSize > 0) {
          vector<char> buffer(bufSize);
          MPI_Recv(&buffer.front(), bufSize, MPI_CHAR, otherProc, 11, Communicator::communicator(), &stat2);
          vector<char>::const_iterator itr = buffer.begin();
          for (unsigned j = 0; j != mSharedNodes[k].size(); ++j) {
            const unsigned i = mSharedNodes[k][j];
            CHECK(i < mNodePositions.size());
            CHECK(nodeOwner.find(i) != nodeOwner.end());
            if (nodeOwner[i] == otherProc) unpackElement(mNodePositions[i], itr, buffer.end());
          }
          CHECK(itr == buffer.end());
        }
      }
    }

    // Wait for our sends to complete.
    if (not sendRequests.empty()) {
      vector<MPI_Status> status(sendRequests.size());
      MPI_Waitall(sendRequests.size(), &sendRequests.front(), &status.front());
    }
  }

  // That's it.
  ENSURE2(validDomainInfo(xmin, xmax, false) == "", validDomainInfo(xmin, xmax, false));
#endif
}
 
//------------------------------------------------------------------------------
// Mesh::generateParallelRind
// Generate a parallel rind of cells around each domain representing a one zone
// thick set of zones shared with the neighboring processors.
// Note we do not recompute the shared elements (nodes & faces) as part of this
// procedure, so following this operation those shared elements are no longer
// on the surface of the local mesh!
//------------------------------------------------------------------------------
// Version without generators.
template<typename Dimension>
void
Mesh<Dimension>::
generateParallelRind() {
  vector<Vector> fakeGenerators(this->numZones());
  vector<SymTensor> fakeHs(this->numZones());
  this->generateParallelRind(fakeGenerators, fakeHs);
}

// Verstion with generators.
template<typename Dimension>
void
Mesh<Dimension>::
generateParallelRind(vector<typename Dimension::Vector>& generators,
                     vector<typename Dimension::SymTensor>& Hs) {

  CONTRACT_VAR(generators);
  CONTRACT_VAR(Hs);

  REQUIRE(generators.size() == this->numZones());
  REQUIRE(Hs.size() == this->numZones());

#ifdef SPHERAL_ENABLE_MPI
  // Parallel procs.
  const unsigned numDomains = Process::getTotalNumberOfProcesses();
  //const unsigned rank = Process::getRank();

  if (numDomains > 1) {
    const unsigned numNeighborDomains = mNeighborDomains.size();

    // Get the bounding coordinates in order to help hashing the coordinates.
    Vector xmin, xmax;
    this->boundingBox(xmin, xmax);

    // Define the hashing scale.
    const double dxhash = (xmax - xmin).maxElement() / double(std::numeric_limits<KeyElement>::max());

    // Puff out the bounds a bit.  We do the all reduce just to ensure
    // bit perfect consistency across processors.
    Vector boxInv;
    for (unsigned i = 0; i != Dimension::nDim; ++i) {
      xmin(i) = allReduce(xmin(i) - dxhash, SPHERAL_OP_MIN);
      xmax(i) = allReduce(xmax(i) + dxhash, SPHERAL_OP_MAX);
      boxInv(i) = safeInv(xmax(i) - xmin(i));
    }

    // Create a lookup for node hashes to IDs.
    typedef boost::bimap<Key, unsigned> Hash2IDType;
    Hash2IDType nodeHash2ID;
    for (unsigned i = 0; i != mNodePositions.size(); ++i) {
      nodeHash2ID.insert(Hash2IDType::value_type(hashPosition(mNodePositions[i], xmin, xmax, boxInv), i));
    }

    // Tell every domain we share a node with about our cells that have that node.
    list<vector<char> > sendBufs;
    vector<unsigned> sendSizes(numNeighborDomains);
    vector<MPI_Request> sendRequests(2*numNeighborDomains);
    for (unsigned kdomain = 0; kdomain != numNeighborDomains; ++kdomain) {
      const unsigned otherProc = mNeighborDomains[kdomain];
      CHECK(mSharedNodes[kdomain].size() > 0);

      // Identify the unique set of cells we need to send.
      vector<unsigned> sendCells;
      for (vector<unsigned>::const_iterator nodeItr = mSharedNodes[kdomain].begin();
           nodeItr != mSharedNodes[kdomain].end();
           ++nodeItr) {
        const vector<unsigned>& cells = this->node(*nodeItr).zoneIDs();
        copy(cells.begin(), cells.end(), back_inserter(sendCells));
      }
      sort(sendCells.begin(), sendCells.end());
      sendCells.erase(unique(sendCells.begin(), sendCells.end()), sendCells.end());

      // Pack up our cells for the other domain.
      sendBufs.push_back(vector<char>());
      vector<char>& buf = sendBufs.back();
      for (vector<unsigned>::const_iterator cellItr = sendCells.begin();
           cellItr != sendCells.end();
           ++cellItr) {
        if (*cellItr != Mesh<Dimension>::UNSETID) {
          packElement(generators[*cellItr], buf);
          packElement(Hs[*cellItr], buf);
          const vector<unsigned>& nodes = mZones[*cellItr].nodeIDs();
          const vector<int>& faces = mZones[*cellItr].faceIDs();
          packElement(unsigned(nodes.size()), buf);
          packElement(unsigned(faces.size()), buf);
          map<unsigned, unsigned> nodeMap;
          for (unsigned inode = 0; inode != nodes.size(); ++inode) {
            packElement(nodeHash2ID.right.at(nodes[inode]), buf);
            nodeMap[nodes[inode]] = inode;
          }
          CHECK(nodeMap.size() == nodes.size());
          for (vector<int>::const_iterator faceItr = faces.begin();
               faceItr != faces.end();
               ++faceItr) {
            const unsigned face = positiveID(*faceItr);
            const vector<unsigned>& faceNodes = mFaces[face].nodeIDs();
            packElement(unsigned(faceNodes.size()), buf);
            if (*faceItr < 0) {
              for (vector<unsigned>::const_reverse_iterator itr = faceNodes.rbegin();
                   itr != faceNodes.rend();
                   ++itr) packElement(nodeMap[*itr], buf);
            } else {
              for (vector<unsigned>::const_iterator itr = faceNodes.begin();
                   itr != faceNodes.end();
                   ++itr) packElement(nodeMap[*itr], buf);
            }
          }
        }
      }
      sendSizes[kdomain] = buf.size();
      MPI_Isend(&sendSizes[kdomain], 1, MPI_UNSIGNED, otherProc, 10, Communicator::communicator(), &sendRequests[2*kdomain]);
      MPI_Isend(&buf.front(), buf.size(), MPI_CHAR, otherProc, 11, Communicator::communicator(), &sendRequests[2*kdomain+1]);
    }
    CHECK(sendBufs.size() == numNeighborDomains);

    // MPI_Barrier(Communicator::communicator());
    // for (unsigned irank = 0; irank != numDomains; ++irank) {
    //   if (rank == irank) {
    //     cerr << "================================================================================" << endl
    //          << "Receiving for domain " << rank << endl;

    // Gather up the cells from our neighboring processors and add them to the local
    // Mesh.
    vector<vector<vector<unsigned> > > newCells;
    for (unsigned kdomain = 0; kdomain != numNeighborDomains; ++kdomain) {
      const unsigned otherProc = mNeighborDomains[kdomain];
      // cerr << "  -----> from domain " << otherProc << endl;
      CHECK(mSharedNodes[kdomain].size() > 0);
      MPI_Status status1, status2;
      unsigned bufSize;
      MPI_Recv(&bufSize, 1, MPI_UNSIGNED, otherProc, 10, Communicator::communicator(), &status1);
      CHECK(bufSize > 0);
      vector<char> buffer(bufSize);
      MPI_Recv(&buffer.front(), bufSize, MPI_CHAR, otherProc, 11, Communicator::communicator(), &status2);
      vector<char>::const_iterator bufItr = buffer.begin();

      // Get the number of nodes and faces for this cell.
      while (bufItr != buffer.end()) {
        generators.push_back(Vector());
        Hs.push_back(SymTensor());
        unsigned numCellNodes, numCellFaces;
        unpackElement(generators.back(), bufItr, buffer.end());
        unpackElement(Hs.back(), bufItr, buffer.end());
        unpackElement(numCellNodes, bufItr, buffer.end());
        unpackElement(numCellFaces, bufItr, buffer.end());

        // Unpack the encoded node positions.
        vector<unsigned> cellNodes;
        cellNodes.reserve(numCellNodes);
        Key hashi;
        for (unsigned i = 0; i != numCellNodes; ++i) {
          unpackElement(hashi, bufItr, buffer.end());
          if (nodeHash2ID.left.find(hashi) == nodeHash2ID.left.end()) {
            nodeHash2ID.insert(Hash2IDType::value_type(hashi, mNodePositions.size()));
            mNodePositions.push_back(quantizedPosition(hashi, xmin, xmax));
          }
          cellNodes.push_back(nodeHash2ID.left.at(hashi));
          // cerr << "   Received node " << cellNodes.back() << " @ " << mNodePositions[cellNodes.back()] << endl;
        }
        CHECK(cellNodes.size() == numCellNodes);
        
        // Unpack the faces for this cell as collections of nodes.
        unsigned nNodesInFace, inode;
        newCells.push_back(vector<vector<unsigned> >());
        for (unsigned k = 0; k != numCellFaces; ++k) {
          unpackElement(nNodesInFace, bufItr, buffer.end());
          newCells.back().push_back(vector<unsigned>(nNodesInFace));
          for (unsigned j = 0; j != nNodesInFace; ++j) {
            unpackElement(inode, bufItr, buffer.end());
            CHECK(inode < cellNodes.size());
            newCells.back().back()[j] = cellNodes[inode];
          }
          CHECK(newCells.back().back().size() == nNodesInFace);
        }
        CHECK(newCells.back().size() == numCellFaces);
      }
    }

    // At this point we have created any necessary new node positions,
    // but have not created the new mesh elements yet.
    // Delegate this stage to the Dimension specific implementations.
    this->createNewMeshElements(newCells);

    //   }
    //   MPI_Barrier(Communicator::communicator());
    // }

    // Wait until all our sends have been satisfied.
    if (not sendRequests.empty()) {
      vector<MPI_Status> status(sendRequests.size());
      MPI_Waitall(sendRequests.size(), &sendRequests.front(), &status.front());
    }
  }
#endif
  ENSURE(generators.size() == this->numZones());
  ENSURE(Hs.size() == this->numZones());
}

//------------------------------------------------------------------------------
// Find unique global IDs for all mesh nodes.  This requires that the 
// Mesh::generateDomainInfo method already be called.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<unsigned> 
Mesh<Dimension>::
globalMeshNodeIDs() const {
  const unsigned numDomains = Process::getTotalNumberOfProcesses();
  const unsigned nlocal = this->numNodes();
  // Pre-conditions.
  VERIFY2((numDomains == 1) or
          (mNeighborDomains.size() > 0 and mSharedNodes.size() == mNeighborDomains.size()),
          "You must call Mesh::generateDomainInfo before calling Mesh::globalMeshNodeIDs: "
          << numDomains << " " << mNeighborDomains.size() << " " << mSharedNodes.size());

  // If we're serial this is easy!
  vector<unsigned> result(nlocal, UNSETID);
#ifdef SPHERAL_ENABLE_MPI
  if (numDomains != 1) {
    unsigned nown;
    const unsigned rank = Process::getRank();
  
    // The parallel case.  Start by having everyone figure out how many nodes
    // they own.
    const unsigned numNeighborDomains = mSharedNodes.size();
    vector<unsigned> owner(nlocal, rank);
    for (unsigned k = 0; k != numNeighborDomains; ++k) {
      CHECK(mNeighborDomains[k] < numDomains);
      if (mNeighborDomains[k] < rank) {
        for (const unsigned i: mSharedNodes[k]) {
          CHECK(i < nlocal);
          owner[i] = mNeighborDomains[k];
        }
      }
    }
    CHECK(owner.size() == nlocal);
    nown = count(owner.begin(), owner.end(), rank);

    // Figure out the range of global IDs for each domain to assign.  This is simple but
    // terribly serial!
    unsigned minID = 0;
    {
      if (rank > 0) {
        MPI_Status recvStatus;
        MPI_Recv(&minID, 1, MPI_UNSIGNED, rank - 1, 10, Communicator::communicator(), &recvStatus);
      }
      if (rank < numDomains - 1) {
        unsigned maxID = minID + nown;
        MPI_Send(&maxID, 1, MPI_UNSIGNED, rank + 1, 10, Communicator::communicator());
      }
    }

    // Assign the global IDs for the vertices we own.
    unsigned j = minID;
    for (unsigned i = 0; i != nlocal; ++i) {
      if (owner[i] == rank) result[i] = j++;
    }
    CHECK(j == minID + nown);

    // Exchange the known IDs between neighboring domains.
    list<vector<unsigned> > sendBuffers, recvBuffers;
    vector<MPI_Request> sendRequests, recvRequests;
    sendRequests.reserve(numNeighborDomains);
    recvRequests.reserve(numNeighborDomains);
    for (unsigned k = 0; k != numNeighborDomains; ++k) {
      const unsigned otherProc = mNeighborDomains[k];
      CHECK(otherProc < numDomains);
      const unsigned n = mSharedNodes[k].size();
      if (mNeighborDomains[k] < rank) {
        recvBuffers.push_back(vector<unsigned>(n, UNSETID));
        recvRequests.push_back(MPI_Request());
        MPI_Irecv(&(recvBuffers.back().front()), n, MPI_UNSIGNED, otherProc, 20, Communicator::communicator(), &recvRequests.back());
      } else {
        sendBuffers.push_back(vector<unsigned>(n, UNSETID));
        sendRequests.push_back(MPI_Request());
        vector<unsigned>& buf = sendBuffers.back();
        for (unsigned i = 0; i != n; ++i) buf[i] = result[mSharedNodes[k][i]];
        MPI_Isend(&buf.front(), n, MPI_UNSIGNED, otherProc, 20, Communicator::communicator(), &sendRequests.back());
      }
    }
    CHECK(sendBuffers.size() == sendRequests.size());
    CHECK(sendRequests.size() <= numNeighborDomains);
    CHECK(recvBuffers.size() == recvRequests.size());
    CHECK(recvRequests.size() <= numNeighborDomains);
    CHECK(sendRequests.size() + recvRequests.size() == numNeighborDomains);

    // Are we receiving anything?
    if (recvRequests.size() > 0) {

      // Wait 'til we have all our receives.
      if (not recvRequests.empty()) {
        vector<MPI_Status> recvStatus(recvRequests.size());
        MPI_Waitall(recvRequests.size(), &recvRequests.front(), &recvStatus.front());
      }

      // Walk the neighbor domains.
      list<vector<unsigned> >::const_iterator recvBufferItr = recvBuffers.begin();
      for (unsigned k = 0; k != numNeighborDomains; ++k) {
        if (mNeighborDomains[k] < rank) {

          // Unpack the global IDs, with the knowledge that a given node can be shared with more
          // than one domain so we have to take a minimum here.
          const unsigned n = mSharedNodes[k].size();
          CHECK(recvBufferItr != recvBuffers.end());
          const vector<unsigned>& buf = *recvBufferItr++;
          CHECK(buf.size() == n);
          for (unsigned j = 0; j != n; ++j) {
            const unsigned i = mSharedNodes[k][j];
            CHECK(i < nlocal);
            result[i] = min(result[i], buf[j]);
          }
        }
      }
      CHECK(recvBufferItr == recvBuffers.end());
    }

    // Don't exit until all of our sends are complete.
    if (not sendRequests.empty()) {
      vector<MPI_Status> sendStatus(sendRequests.size());
      MPI_Waitall(sendRequests.size(), &sendRequests.front(), &sendStatus.front());
    }

    // I think this post-condition check is too expensive to always have in a 
    // contract, so I'm putting this same logic in the Mesh unit tests instead.
//     // Post-conditions (parallel).
//     BEGIN_CONTRACT_SCOPE
//     {
//       // Make sure everyone is consistent about the global IDs of shared nodes.
//       list<vector<unsigned> > sendBuffers, recvBuffers;
//       vector<MPI_Request> sendRequests, recvRequests;
//       sendRequests.reserve(numNeighborDomains);
//       recvRequests.reserve(numNeighborDomains);
//       for (unsigned k = 0; k != numNeighborDomains; ++k) {
//         const unsigned otherProc = mNeighborDomains[k];
//         CHECK(otherProc < numDomains);
//         const unsigned n = mSharedNodes[k].size();

//         recvBuffers.push_back(vector<unsigned>(n, UNSETID));
//         recvRequests.push_back(MPI_Request());
//         MPI_Irecv(&(recvBuffers.back().front()), n, MPI_UNSIGNED, otherProc, 20, Communicator::communicator(), &recvRequests.back());

//         sendBuffers.push_back(vector<unsigned>(n, UNSETID));
//         sendRequests.push_back(MPI_Request());
//         vector<unsigned>& buf = sendBuffers.back();
//         for (unsigned i = 0; i != n; ++i) buf[i] = result[mSharedNodes[k][i]];
//         MPI_Isend(&buf.front(), n, MPI_UNSIGNED, otherProc, 20, Communicator::communicator(), &sendRequests.back());
//       }
//       CHECK(sendBuffers.size() == numNeighborDomains);
//       CHECK(sendRequests.size() == numNeighborDomains);
//       CHECK(recvBuffers.size() == numNeighborDomains);
//       CHECK(recvRequests.size() == numNeighborDomains);

//       // Wait 'til we have all our receives.
//       vector<MPI_Status> recvStatus(recvRequests.size());
//       MPI_Waitall(recvRequests.size(), &recvRequests.front(), &recvStatus.front());

//       // Walk the neighbor domains.
//       list<vector<unsigned> >::const_iterator recvBufferItr = recvBuffers.begin();
//       for (unsigned k = 0; k != numNeighborDomains; ++k) {
//         const unsigned n = mSharedNodes[k].size();
//         CHECK(recvBufferItr != recvBuffers.end());
//         const vector<unsigned>& buf = *recvBufferItr++;
//         CHECK(buf.size() == n);
//         for (unsigned j = 0; j != n; ++j) {
//           const unsigned i = mSharedNodes[k][j];
//           CHECK(i < nlocal);
//           ENSURE(result[i] == buf[j]);
//         }
//       }
//       CHECK(recvBufferItr == recvBuffers.end());

//       // Don't exit until all of our sends are complete.
//       if (sendRequests.size() > 0) {
//         vector<MPI_Status> sendStatus(sendRequests.size());
//         MPI_Waitall(sendRequests.size(), &sendRequests.front(), &sendStatus.front());
//       }
//     }
//     END_CONTRACT_SCOPE
  }
#endif


  // Post-conditions.
  ENSURE(result.size() == nlocal);
  ENSURE(count(result.begin(), result.end(), UNSETID) == 0);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Define unique global IDs for all mesh faces.
// Note unlike the node IDs the result of this method is not numbered 
// sequentially across domains.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<unsigned> 
Mesh<Dimension>::
globalMeshFaceIDs(const vector<unsigned>& globalNodeIDs) const {

  // Now hash the faceIDs based on the global nodes. 
  vector<unsigned> result;
  result.reserve(mFaces.size());
  for (const Face& face: mFaces) {
    const vector<unsigned>& locals = face.nodeIDs();
    vector<unsigned> globals;
    globals.reserve(locals.size());
    for (const unsigned i: locals) globals.push_back(globalNodeIDs[i]);
    sort(globals.begin(), globals.end());
    size_t seed = 0;
    for (const unsigned i: globals) boost::hash_combine(seed, i);
    result.push_back(seed);
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    ENSURE(result.size() == mFaces.size());

    // Make sure the local results are unique.
    for (int i = 0; i < (int)result.size() - 1; ++i) {
      ENSURE(find(result.begin() + i + 1, result.end(), result[i]) == result.end());
    }
  }
  END_CONTRACT_SCOPE

  return result;
}

//------------------------------------------------------------------------------
// Fill in the NodeList offsets.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Mesh<Dimension>::
storeNodeListOffsets(const vector<NodeList<Dimension>*>& nodeListPtrs,
                     const vector<unsigned>& offsets) {
  VERIFY(nodeListPtrs.size() == offsets.size());
  this->storeNodeListOffsets(nodeListPtrs.begin(), nodeListPtrs.end(), offsets);
}

//------------------------------------------------------------------------------
// Find the minimum bounding box for the mesh.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Mesh<Dimension>::
boundingBox(typename Dimension::Vector& xmin,
            typename Dimension::Vector& xmax) const {
  Spheral::boundingBox(mNodePositions, xmin, xmax);
#ifdef SPHERAL_ENABLE_MPI
  for (unsigned i = 0; i != Dimension::nDim; ++i) {
    xmin(i) = allReduce(xmin(i), SPHERAL_OP_MIN);
    xmax(i) = allReduce(xmax(i), SPHERAL_OP_MAX);
  }
#endif
}

// //------------------------------------------------------------------------------
// // Internal method to fill in extra comm data based on the shared node info.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void
// Mesh<Dimension>::
// buildAncillaryCommData() {

//   // Clear out old data.
//   cerr << Process::getRank() << " Mesh::buildAncillaryCommData 1" << endl;
//   mSharedFaces = vector<vector<unsigned> >();

// #ifdef SPHERAL_ENABLE_MPI
//   // Compute the shared faces.  First get the unique global face IDs.
//   const vector<unsigned> globalNodeIDs = this->globalMeshNodeIDs();
//   const vector<unsigned> globalFaceIDs = this->globalMeshFaceIDs(globalNodeIDs);

//   // Compute the inverse face id mapping: global->local.
//   cerr << Process::getRank() << " Mesh::buildAncillaryCommData 2" << endl;
//   map<unsigned, unsigned> global2local;
//   for (unsigned i = 0; i != globalFaceIDs.size(); ++i) global2local[globalFaceIDs[i]] = i;
//   CHECK(globalFaceIDs.size() == global2local.size());

//   // Make a sorted version of the global face IDs.
//   cerr << Process::getRank() << " Mesh::buildAncillaryCommData 3" << endl;
//   vector<unsigned> sortedGlobalIDs(globalFaceIDs);
//   sort(sortedGlobalIDs.begin(), sortedGlobalIDs.end());

//   // Send our global face IDs to each of our neighbors.
//   unsigned numFaces = globalFaceIDs.size();
//   vector<MPI_Request> requests;
//   requests.reserve(2*mNeighborDomains.size());
//   for (unsigned neighborProc: mNeighborDomains) {
//     requests.push_back(MPI_Request());
//     MPI_Isend(&numFaces, 1, MPI_UNSIGNED, neighborProc, 1, Communicator::communicator(), &requests.back());
//     if (numFaces > 0) {
//       requests.push_back(MPI_Request());
//       MPI_Isend(&sortedGlobalIDs.front(), numFaces, MPI_UNSIGNED, neighborProc, 2, Communicator::communicator(), &requests.back());
//     }
//   }
//   CHECK(requests.size() <= 2*mNeighborDomains.size());

//   // Now go through each of our neighbors shared faces, and build the set we share with 
//   // them.  We store the local face indices in the order of the sorted global face IDs,
//   // which should ensure that both processors agree on the face order.
//   cerr << Process::getRank() << " Mesh::buildAncillaryCommData 4" << endl;
//   for (unsigned neighborProc: mNeighborDomains) {
//     unsigned numOtherFaces;
//     MPI_Status recvStatus;
//     MPI_Recv(&numOtherFaces, 1, MPI_UNSIGNED, neighborProc, 1, Communicator::communicator(), &recvStatus);
//     mSharedFaces.push_back(vector<unsigned>());
//     if (numOtherFaces > 0) {
//       vector<unsigned> otherGlobalFaceIDs(numOtherFaces);
//       MPI_Recv(&otherGlobalFaceIDs.front(), numOtherFaces, MPI_UNSIGNED, neighborProc, 2, Communicator::communicator(), &recvStatus);
//       vector<unsigned>::iterator lastItr = sortedGlobalIDs.begin(), itr;
//       for (unsigned iglobal: otherGlobalFaceIDs) {
//         itr = lower_bound(lastItr, sortedGlobalIDs.end(), iglobal);
//         if (itr != sortedGlobalIDs.end()) {
//           mSharedFaces.back().push_back(global2local[iglobal]);
//           lastItr = itr;
//         }
//       }
//     }
//   }
//   CHECK(mSharedFaces.size() == mNeighborDomains.size());

//   // Wait until all our sends are completed.
//   vector<MPI_Status> recvStatus(requests.size());
//   MPI_Waitall(requests.size(), &requests.front(), &recvStatus.front());
//   cerr << Process::getRank() << " Mesh::buildAncillaryCommData 5" << endl;
// #endif

//   // Post-conditions.
//   ENSURE(mSharedNodes.size() == mNeighborDomains.size());
//   ENSURE(mSharedFaces.size() == mNeighborDomains.size());
// #ifdef SPHERAL_ENABLE_MPI
//   BEGIN_CONTRACT_SCOPE
//   {
//     cerr << Process::getRank() << " Mesh::buildAncillaryCommData 6" << endl;
//     // Make sure each domain agrees about the shared faces.
//     vector<MPI_Request> requests;
//     requests.reserve(2*mNeighborDomains.size());
//     vector<unsigned> numSharedFaces(mNeighborDomains.size());
//     vector<vector<unsigned> > globalSharedFaces(mNeighborDomains.size());
//     for (unsigned idomain = 0; idomain != mNeighborDomains.size(); ++idomain) {
//       for (unsigned i: mSharedFaces[idomain]) globalSharedFaces[idomain].push_back(globalFaceIDs[i]);
//     }

//     cerr << Process::getRank() << "Mesh::buildAncillaryCommData 7" << endl;
//     for (unsigned idomain = 0; idomain != mNeighborDomains.size(); ++idomain) {
//       requests.push_back(MPI_Request());
//       numSharedFaces[idomain] = mSharedFaces[idomain].size();
//       CHECK(globalSharedFaces[idomain].size() == numSharedFaces[idomain]);
//       MPI_Isend(&numSharedFaces[idomain], 1, MPI_UNSIGNED, mNeighborDomains[idomain], 1, Communicator::communicator(), &requests.back());
//       if (numSharedFaces[idomain] > 0) {
//         requests.push_back(MPI_Request());
//         MPI_Isend(&globalSharedFaces[idomain].front(), numSharedFaces[idomain], MPI_UNSIGNED, mNeighborDomains[idomain], 2, Communicator::communicator(), &requests.back());
//       }
//     }
//     CHECK(requests.size() <= 2*mNeighborDomains.size());

//     cerr << Process::getRank() << " Mesh::buildAncillaryCommData 8" << endl;
//     for (unsigned idomain = 0; idomain != mNeighborDomains.size(); ++idomain) {
//       unsigned numOtherFaces;
//       MPI_Status recvStatus;
//       MPI_Recv(&numOtherFaces, 1, MPI_UNSIGNED, mNeighborDomains[idomain], 1, Communicator::communicator(), &recvStatus);
//       ENSURE(numOtherFaces == mSharedFaces[idomain].size());
//       if (numOtherFaces > 0) {
//         vector<unsigned> otherFaces(numOtherFaces);
//         MPI_Recv(&otherFaces.front(), numOtherFaces, MPI_UNSIGNED, mNeighborDomains[idomain], 2, Communicator::communicator(), &recvStatus);
//         ENSURE(otherFaces == globalSharedFaces[idomain]);
//       }
//     }

//     vector<MPI_Status> recvStats(requests.size());
//     MPI_Waitall(requests.size(), &requests.front(), &recvStatus.front());
//   }
//   END_CONTRACT_SCOPE
// #endif
//   cerr << Process::getRank() << " Mesh::buildAncillaryCommData 9" << endl;
// }

//------------------------------------------------------------------------------
// Basic mesh validity checks.
//------------------------------------------------------------------------------
template<typename Dimension>
string
Mesh<Dimension>::
valid() const {
  std::stringstream result;

  // Check that each face has two cells, and has orientation.
  for (const Face& face: mFaces) {
    if (not ((face.zone1ID() <  0 and face.zone2ID() >= 0) or
             (face.zone1ID() >= 0 and face.zone2ID() <  0))) {
      result << "Expected one negative zone ID for face " 
             << face.ID() << " : " << face.zone1ID() << " " << face.zone2ID();
      return result.str();
    }
  }

  return result.str();
}
//------------------------------------------------------------------------------
// Check that the communicated information is consistent.
// Optionally we can also check that each send node is owned by only one 
// process.
//------------------------------------------------------------------------------
template<typename Dimension>
string
Mesh<Dimension>::
validDomainInfo(const typename Dimension::Vector& xmin,
                const typename Dimension::Vector& xmax,
                const bool /*checkUniqueSendProc*/) const {
  string result = "";

  CONTRACT_VAR(xmin);
  CONTRACT_VAR(xmax);
#ifdef SPHERAL_ENABLE_MPI
  const unsigned rank = Process::getRank();
  const unsigned numDomains = Process::getTotalNumberOfProcesses();

  // The inverse box scale.
  Vector boxInv;
  for (unsigned i = 0; i != Dimension::nDim; ++i) boxInv(i) = safeInv(xmax(i) - xmin(i));

  // First check that all processors agree about who is talking to whom.
  {
    for (unsigned sendProc = 0; sendProc != numDomains; ++sendProc) {
      unsigned numNeighbors = mNeighborDomains.size();
      vector<unsigned> checkNeighbors = mNeighborDomains;
      MPI_Bcast(&numNeighbors, 1, MPI_UNSIGNED, sendProc, Communicator::communicator());
      if (numNeighbors > 0) {
        checkNeighbors.resize(numNeighbors);
        MPI_Bcast(&checkNeighbors.front(), numNeighbors, MPI_UNSIGNED, sendProc, Communicator::communicator());
        if (not (binary_search(checkNeighbors.begin(), checkNeighbors.end(), rank) == 
                 binary_search(mNeighborDomains.begin(), mNeighborDomains.end(), sendProc))) {
          result = "Processors don't agree about who is talking to whom!";
        }
      }
    }
    result = reduceToMaxString(result, rank, numDomains);
  }

  // Check that the processors that are talking to agree about the number of shared nodes.
  if (result == "") {
    CHECK(mSharedNodes.size() == mNeighborDomains.size());
    vector<MPI_Request> requests(2*mNeighborDomains.size());
    vector<unsigned> numLocalSharedNodes, numOtherSharedNodes(mNeighborDomains.size());
    for (unsigned k = 0; k != mSharedNodes.size(); ++k) numLocalSharedNodes.push_back(mSharedNodes[k].size());
    CHECK(numLocalSharedNodes.size() == mNeighborDomains.size());
    CHECK(numOtherSharedNodes.size() == mNeighborDomains.size());
    for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
      const unsigned otherProc = mNeighborDomains[k];
      MPI_Isend(&numLocalSharedNodes[k], 1, MPI_UNSIGNED, otherProc,
                (rank + 1)*numDomains + otherProc, Communicator::communicator(), &requests[k]);
      MPI_Irecv(&numOtherSharedNodes[k], 1, MPI_UNSIGNED, otherProc,
                (otherProc + 1)*numDomains + rank, Communicator::communicator(), &requests[mNeighborDomains.size() + k]);
    }
    if (not requests.empty()) {
      vector<MPI_Status> status(requests.size());
      MPI_Waitall(requests.size(), &requests.front(), &status.front());
    }
    if (not (numLocalSharedNodes == numOtherSharedNodes)) {
      result = "Processors don't agree about number of shared nodes.";
    }
    result = reduceToMaxString(result, rank, numDomains);
  }
      
  // Similary check the number of shared faces.
  if (result == "") {
    CHECK(mSharedFaces.size() == mNeighborDomains.size());
    vector<MPI_Request> requests(2*mNeighborDomains.size());
    vector<unsigned> numLocalSharedFaces, numOtherSharedFaces(mNeighborDomains.size());
    for (unsigned k = 0; k != mSharedFaces.size(); ++k) numLocalSharedFaces.push_back(mSharedFaces[k].size());
    CHECK(numLocalSharedFaces.size() == mNeighborDomains.size());
    CHECK(numOtherSharedFaces.size() == mNeighborDomains.size());
    for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
      const unsigned otherProc = mNeighborDomains[k];
      MPI_Isend(&numLocalSharedFaces[k], 1, MPI_UNSIGNED, otherProc,
                (rank + 1)*numDomains + otherProc, Communicator::communicator(), &requests[k]);
      MPI_Irecv(&numOtherSharedFaces[k], 1, MPI_UNSIGNED, otherProc,
                (otherProc + 1)*numDomains + rank, Communicator::communicator(), &requests[mNeighborDomains.size() + k]);
    }
    if (not requests.empty()) {
      vector<MPI_Status> status(requests.size());
      MPI_Waitall(requests.size(), &requests.front(), &status.front());
    }
    if (not (numLocalSharedFaces == numOtherSharedFaces)) {
      result = "Processors don't agree about number of shared faces.";
    }
    result = reduceToMaxString(result, rank, numDomains);
  }
      
  // For now we are suspending the hashing checks due to some kind of accuracy conflict with polytope.

  // // Check that we agree with all our neighbors about which nodes we share.
  // if (result == "") {
  //   vector<Key> allLocalHashes;
  //   vector<vector<Key> > otherHashes;
  //   for (unsigned i = 0; i != mNodePositions.size(); ++i) 
  //     allLocalHashes.push_back(hashPosition(mNodePositions[i], xmin, xmax, boxInv));
  //   CHECK(allLocalHashes.size() == mNodePositions.size());
  //   exchangeTuples(allLocalHashes, mNeighborDomains, mSharedNodes, otherHashes);

  //   for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
  //     for (unsigned i = 0; i != mSharedNodes[k].size(); ++i) {
  //       if (not (allLocalHashes[mSharedNodes[k][i]] == otherHashes[k][i]))
  //         result = "Node hash neighbor comparisons fail!";
  //     }
  //   }
  //   result = reduceToMaxString(result, rank, numDomains);
  // }

  // // Ditto for shared faces.
  // if (result == "") {
  //   vector<Key> allLocalHashes;
  //   vector<vector<Key> > otherHashes;
  //   for (unsigned i = 0; i != mFaces.size(); ++i) 
  //     allLocalHashes.push_back(hashPosition(mFaces[i].position(), xmin, xmax, boxInv));
  //   CHECK(allLocalHashes.size() == mFaces.size());
  //   exchangeTuples(allLocalHashes, mNeighborDomains, mSharedFaces, otherHashes);

  //   for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
  //     for (unsigned i = 0; i != mSharedFaces[k].size(); ++i) {
  //       if (not (allLocalHashes[mSharedFaces[k][i]] == otherHashes[k][i]))
  //         result = "Face hash neighbor comparisons fail!";
  //     }
  //   }
  //   result = reduceToMaxString(result, rank, numDomains);
  // }

  // // Check the uniqueness of shared nodes (only shared with one other domain unless owned).
  // if (result == "" and checkUniqueSendProc) {
  //   unordered_map<unsigned, unsigned> shareCount;
  //   unordered_map<unsigned, vector<unsigned> > nodeSharedWithDomains;
  //   for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
  //     const unsigned otherProc = mNeighborDomains[k];
  //     CHECK(k < mSharedNodes.size());
  //     for (unsigned j = 0; j != mSharedNodes[k].size(); ++j) {
  //       CHECK(j < mSharedNodes[k].size());
  //       const unsigned i = mSharedNodes[k][j];
  //       nodeSharedWithDomains[i].push_back(otherProc);
  //     }
  //   }
  //   for (unordered_map<unsigned, vector<unsigned> >::const_iterator itr = nodeSharedWithDomains.begin();
  //        itr != nodeSharedWithDomains.end();
  //        ++itr) {
  //     const unsigned i = itr->first;
  //     const vector<unsigned>& otherDomains = itr->second;
  //     if (otherDomains.size() > 1 and
  //         *min_element(otherDomains.begin(), otherDomains.end()) < rank) {
  //       stringstream thpt;
  //       thpt << " Node " << i << " on domain " << rank << " shared with domains [";
  //       for (unsigned k = 0; k != otherDomains.size(); ++k) thpt << " " << otherDomains[k];
  //       thpt << "]" << endl;
  //       result += thpt.str();
  //     }
  //   }
  //   result = reduceToMaxString(result, rank, numDomains);
  // }
#endif

  return result;
}

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<typename Dimension>
const unsigned Mesh<Dimension>::UNSETID = std::numeric_limits<int>::max();

}
