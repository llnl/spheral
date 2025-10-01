//---------------------------------Spheral++----------------------------------//
// LineMesh -- 1-D mesh class.
//
// Created by JMO, Tue Oct 12 23:07:22 PDT 2010
//----------------------------------------------------------------------------//
#include "Mesh.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/bisectSearch.hh"
#include "Distributed/allReduce.hh"
#include "Distributed/Communicator.hh"
#include "Utilities/DBC.hh"

#include "config.hh"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <limits>
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
// Comparator function to sort indices.
//------------------------------------------------------------------------------
struct CompareIndicesByPositions {
  CompareIndicesByPositions(const vector<Dim<1>::Vector>& positions):
    mPositions(positions) {}
  bool operator()(const unsigned lhs, const unsigned rhs) {
    REQUIRE(lhs < mPositions.size());
    REQUIRE(rhs < mPositions.size());
    return mPositions[lhs].x() < mPositions[rhs].x();
  }
  const vector<Dim<1>::Vector>& mPositions;
};

//------------------------------------------------------------------------------
// Mesh::reconstruct(generators, xmin, xmax)
//------------------------------------------------------------------------------
template<>
void
Mesh<Dim<1> >::
reconstructInternal(const vector<Mesh<Dim<1> >::Vector>& localGenerators,
                    const Mesh<Dim<1> >::Vector& xmin,
                    const Mesh<Dim<1> >::Vector& xmax) {

  // Is there anything to do?
  if (allReduce(unsigned(localGenerators.size()), SPHERAL_OP_SUM) == 0) return;

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(xmin < xmax);
    for (vector<Vector>::const_iterator itr = localGenerators.begin();
         itr != localGenerators.end();
         ++itr) {
      REQUIRE2(xmin <= *itr and *itr <= xmax, "Node out of bounds:  " << *itr << " not in [" << xmin << " " << xmax << "]");
    }
  }
  END_CONTRACT_SCOPE

  // Get the full set of generators we need.
  vector<Vector> generators = localGenerators;
#ifdef ENABLE_MPI
  // Parallel info.
  const unsigned rank = Process::getRank();
  const unsigned numDomains = Process::getTotalNumberOfProcesses();
  // We need the parallel sets of generators.
  {  
    // Create the set of hashed local generators.
    // const double degeneracy = 1.0e-8*max(1.0, xmax.x() - xmin.x());
    // const double dxInv = 1.0/degeneracy;
    // vector<KeyElement> localHashes;
    // localHashes.reserve(localGenerators.size());
    // for (unsigned i = 0; i != localGenerators.size(); ++i) 
    //   localHashes.push_back(KeyElement(dxInv*(localGenerators[i].x() - xmin.x()) + 0.5));

    const Vector boxInv(safeInv(xmax.x() - xmin.x()));
    set<Key> occupiedHashes;
    for (unsigned i = 0; i != localGenerators.size(); ++i) 
      occupiedHashes.insert(hashPosition(localGenerators[i], xmin, xmax, boxInv));
    CHECK(occupiedHashes.size() == localGenerators.size());

    vector<char> localBuffer;
    Key hashi;
    Vector xi;
    for (vector<Vector>::const_iterator itr = localGenerators.begin();
         itr != localGenerators.end();
         ++itr) packElement(*itr, localBuffer);
    for (unsigned sendProc = 0; sendProc != numDomains; ++sendProc) {
      vector<char> buffer = localBuffer;
      unsigned bufSize = buffer.size();
      MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, Communicator::communicator());
      if (bufSize > 0) {
        buffer.resize(bufSize);
        MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, Communicator::communicator());
        if (sendProc != rank) {
          vector<char>::const_iterator bufItr = buffer.begin();
          while (bufItr != buffer.end()) {
            unpackElement(xi, bufItr, buffer.end());
            hashi = hashPosition(xi, xmin, xmax, boxInv);
            if (occupiedHashes.find(hashi) == occupiedHashes.end()) {
              generators.push_back(xi);
              occupiedHashes.insert(hashi);
            }
          }
        }
      }
    }
  }
#endif

  // Find the sorted order for the zones by index.
  vector<unsigned> zoneOrder;
  for (unsigned i = 0; i != generators.size(); ++i) zoneOrder.push_back(i);
  CompareIndicesByPositions zoneComparator(generators);
  sort(zoneOrder.begin(), zoneOrder.end(), zoneComparator);
  CHECK(zoneOrder.size() == generators.size());
  BEGIN_CONTRACT_SCOPE
  {
    // This check not only checks that the positions are sorted but
    // also that there are no duplicates.
    for (vector<unsigned>::const_iterator itr = zoneOrder.begin();
         itr < zoneOrder.end() - 1;
         ++itr) {
      CHECK2(generators[*itr].x() < generators[*(itr + 1)].x(), 
             "Bad generator positions:  " << generators[*itr] << " " << generators[*(itr + 1)]);
    }
  }
  END_CONTRACT_SCOPE

  // It's also useful to have the reverse mapping.
  vector<unsigned> sortedZoneIDs(zoneOrder.size(), UNSETID);
  for (unsigned i = 0; i != zoneOrder.size(); ++i) sortedZoneIDs[zoneOrder[i]] = i;
  CHECK(find(sortedZoneIDs.begin(), sortedZoneIDs.end(), UNSETID) == sortedZoneIDs.end());

  // Compute the node positions from the zones.
  mNodePositions.push_back(xmin);
  for (vector<unsigned>::const_iterator itr = zoneOrder.begin();
       itr < zoneOrder.end() - 1;
       ++itr) mNodePositions.push_back(0.5*(generators[*itr] + generators[*(itr + 1)]));
  mNodePositions.push_back(xmax);
  // cerr << "Sorted generators:  ";
  // for (unsigned i = 0; i != zoneOrder.size(); ++i) cerr << " " << generators[zoneOrder[i]];
  // cerr << endl << "Node positions : ";
  // for (unsigned i = 0; i != mNodePositions.size(); ++i) cerr << " " << mNodePositions[i];
  // cerr << endl;
  CHECK2(mNodePositions.size() == generators.size() + 1, "Blago!  " << mNodePositions.size() << " " << generators.size());

  // Construct the nodes from the positions.
  mNodes.reserve(mNodePositions.size());
  for (unsigned i = 0; i != mNodePositions.size(); ++i) {
    vector<unsigned> zones;
    if (i > 0)                         zones.push_back(zoneOrder[i - 1]);
    if (i < mNodePositions.size() - 1) zones.push_back(zoneOrder[i]);
    mNodes.push_back(Node(*this, i, zones));
    // cerr << "Node -> zones:  " << i << " " << zones[0] << " " << zones[1] << " : " << mNodePositions[i] << " "
    //      << (zones[0] != UNSETID ? generators[zones[0]] : -1.0) << " "
    //      << (zones[1] != UNSETID ? generators[zones[1]] : -1.0) << " "
    //      << endl;
  }
  CHECK(mNodes.size() == mNodePositions.size());

  // Construct the edges.  For a 1-D LineMesh edges are degenerate -- they only
  // have one node.
  mEdges.reserve(mNodePositions.size());
  for (unsigned i = 0; i != mNodePositions.size(); ++i) mEdges.push_back(Edge(*this, i, i, i));
  CHECK(mEdges.size() == mNodePositions.size());

  // Construct the faces.  Same as the edges -- one per node.
  mFaces.reserve(mNodePositions.size());
  mFaces.push_back(Face(*this,
                        0,
                        UNSETID,
                        ~zoneOrder[0],
                        vector<unsigned>(1, 0)));
  for (unsigned i = 1; i < mNodePositions.size() - 1; ++i) mFaces.push_back(Face(*this,
                                                                                 i,
                                                                                 zoneOrder[i - 1],
                                                                                 ~zoneOrder[i],
                                                                                 vector<unsigned>(1, i)));
  mFaces.push_back(Face(*this,
                        mNodePositions.size() - 1,
                        zoneOrder[mNodePositions.size() - 2],
                        ~UNSETID,
                        vector<unsigned>(1, mNodePositions.size() - 1)));
  CHECK(mFaces.size() == mNodePositions.size());
  // for (unsigned i = 0; i != mFaces.size(); ++i) {
  //   cerr << "Face:  " << i << " " << mFaces[i].zone1ID() << " " << mFaces[i].zone2ID() << " : "
  //        << mFaces[i].position() << " " 
  //        << (mFaces[i].zone1ID() < generators.size() ? generators[mFaces[i].zone1ID()] : -1.0) << " "
  //        << (mFaces[i].zone2ID() < generators.size() ? generators[mFaces[i].zone2ID()] : -1.0) << " "
  //        << endl;
  // }

  // Finally construct the zones.
  for (unsigned igen = 0; igen != generators.size(); ++igen) {
    const int node1 = sortedZoneIDs[igen];
    const int node2 = node1 + 1;
    CHECK2(generators[igen].x() > mNodePositions[node1].x() and
           generators[igen].x() < mNodePositions[node2].x(),
           "Generator outside node boundaries!  "
           << igen << " " << node1 << " " << node2 << " : "
           << generators[igen] << " "
           << mNodePositions[node1] << " "
           << mNodePositions[node2]);

    // Now we can build the zone.
    // We use the fact that there is a one to one mapping of nodes->faces.
    vector<int> faceIDs;
    faceIDs.push_back(~node1);
    faceIDs.push_back(node2);
    mZones.push_back(Zone(*this, igen, faceIDs));
  }

  // Remove the zones for generators that are not local.
  vector<unsigned> mask(generators.size(), 0);
  fill(mask.begin(), mask.begin() + localGenerators.size(), 1);
  this->removeZonesByMask(mask);

  // Build the parallel info.
  this->generateDomainInfo();
  // MPI_Barrier(Communicator::communicator());
  // for (unsigned irank = 0; irank != Process::getTotalNumberOfProcesses(); ++irank) {
  //   if (Process::getRank() == irank) {
  //     cerr << "LineMesh neighborDomains : " << mNeighborDomains.size() << " : ";
  //     copy(mNeighborDomains.begin(), mNeighborDomains.end(), ostream_iterator<unsigned>(cerr, " "));
  //     cerr << endl;
  //   }
  //   MPI_Barrier(Communicator::communicator());
  // }

  // Post-conditions.
  ENSURE(mNodePositions.size() == localGenerators.size() + 1);
  ENSURE(mNodes.size() == localGenerators.size() + 1);
  ENSURE(mEdges.size() == localGenerators.size() + 1);
  ENSURE(mFaces.size() == localGenerators.size() + 1);
  ENSURE(mZones.size() == localGenerators.size());
  ENSURE(mSharedNodes.size() == mNeighborDomains.size());
  ENSURE(mSharedFaces.size() == mNeighborDomains.size());
  BEGIN_CONTRACT_SCOPE
  {
    // In 1D we know that every shared node should correspond to a shared face!
    for (unsigned i = 0; i != mSharedNodes.size(); ++i) ENSURE(mSharedFaces[i].size() == mSharedNodes[i].size());
  }
  END_CONTRACT_SCOPE
  ENSURE2(this->valid() == "", this->valid());
}

//------------------------------------------------------------------------------
// Mesh::reconstruct(generators, polygon)
//------------------------------------------------------------------------------
template<>
void
Mesh<Dim<1> >::
reconstruct(const vector<Dim<1>::Vector>& generators,
            const Dim<1>::FacetedVolume& boundary) {
  this->reconstruct(generators, boundary.xmin(), boundary.xmax());
}

//------------------------------------------------------------------------------
// Compute the bounding surface of the mesh.
//------------------------------------------------------------------------------
template<>
Dim<1>::FacetedVolume
Mesh<Dim<1> >::
boundingSurface() const {
  // The best we can do is the bounding vertex positions.
  double xmin = std::numeric_limits<double>::max(), 
         xmax = std::numeric_limits<double>::lowest();
  for (unsigned i = 0; i != mNodePositions.size(); ++i) {
    xmin = std::min(xmin, mNodePositions[i].x());
    xmax = std::max(xmax, mNodePositions[i].x());
  }
  xmin = allReduce(xmin, SPHERAL_OP_MIN);
  xmax = allReduce(xmax, SPHERAL_OP_MAX);
  return FacetedVolume(Vector(0.5*(xmin + xmax)), 0.5*(xmax - xmin));
}

//------------------------------------------------------------------------------
// Add new mesh elements.
//------------------------------------------------------------------------------
template<>
void
Mesh<Dim<1> >::
createNewMeshElements(const vector<vector<vector<unsigned> > >& newCells) {

  // Pre-conditions.
  REQUIRE(mNodes.size() <= mNodePositions.size());
  BEGIN_CONTRACT_SCOPE
  {
    for (const vector<vector<unsigned> >& cellFaces: newCells) {
      CONTRACT_VAR(cellFaces);
      REQUIRE(cellFaces.size() == 2);
      REQUIRE(cellFaces[0].size() == 1 and cellFaces[0][0] < mNodePositions.size());
      REQUIRE(cellFaces[1].size() == 1 and cellFaces[1][0] < mNodePositions.size());
    }
  }
  END_CONTRACT_SCOPE

  // Some useful sizes.
  const unsigned numOldNodes = mNodes.size();
  const unsigned numNewNodes = mNodePositions.size();
  const unsigned numOldZones = mZones.size();

  // Copy the starting connectivity from nodes->zones.
  map<unsigned, set<unsigned> > nodeCells;
  for (unsigned inode = 0; inode != numOldNodes; ++inode) {
    nodeCells[inode] = set<unsigned>(mNodes[inode].mZoneIDs.begin(),
                                     mNodes[inode].mZoneIDs.end());
  }

  // Fill in the existing cell centroids.
  vector<double> cellX(numOldZones + newCells.size());
  for (unsigned i = 0; i != numOldZones; ++i) {
    CHECK(mZones[i].mNodeIDs.size() == 2);
    cellX[i] = 0.5*(mNodePositions[mZones[i].mNodeIDs[0]].x() + mNodePositions[mZones[i].mNodeIDs[1]].x());
  }

  // Update the map of nodes->zones elements with our new cells.
  // We also store the positions of the new cell centroids.
  for (unsigned k = 0; k != newCells.size(); ++k) {
    const unsigned inode1 = newCells[k][0][0];
    const unsigned inode2 = newCells[k][1][0];
    CHECK(inode1 < numNewNodes and
          inode2 < numNewNodes);
    nodeCells[inode1].insert(numOldZones + k);
    nodeCells[inode2].insert(numOldZones + k);
    if (inode1 < numOldNodes) {
      set<unsigned>::iterator itr1 = nodeCells[inode1].find(UNSETID);
      if (itr1 != nodeCells[inode1].end()) nodeCells[inode1].erase(itr1);
    }
    if (inode2 < numOldNodes) {
      set<unsigned>::iterator itr2 = nodeCells[inode2].find(UNSETID);
      if (itr2 != nodeCells[inode2].end()) nodeCells[inode2].erase(itr2);
    }
    cellX[numOldZones + k] = 0.5*(mNodePositions[inode1].x() + mNodePositions[inode2].x());
  }

  // Create the new nodes, edges, and faces.
  for (unsigned inode = numOldNodes; inode != numNewNodes; ++inode) {
    CHECK2(nodeCells[inode].size() == 1 or nodeCells[inode].size() == 2, nodeCells[inode].size() << " : " << inode << " " << mNodes.size() << " " << mNodePositions.size());
    if (nodeCells[inode].size() == 1) nodeCells[inode].insert(UNSETID);
    mNodes.push_back(Node(*this, inode, vector<unsigned>(nodeCells[inode].begin(), nodeCells[inode].end())));
    mEdges.push_back(Edge(*this, inode, inode, inode));

    // For the face we have to figure out which cell is below and above us in x.
    int z1 = min(mNodes[inode].mZoneIDs[0], mNodes[inode].mZoneIDs[1]);
    int z2 = max(mNodes[inode].mZoneIDs[0], mNodes[inode].mZoneIDs[1]);
    CHECK(z1 < z2 and 
          z1 < (int)numOldZones + (int)newCells.size() and 
          (z2 == (int)UNSETID or z2 < (int)numOldZones + (int)newCells.size()));
    if (cellX[z1] > mNodePositions[inode].x()) std::swap(z1, z2);
    mFaces.push_back(Face(*this, inode, ~z1, z2, vector<unsigned>(1, inode)));
  }

  // Create the new cells.
  for (unsigned k = 0; k != newCells.size(); ++k) {
    vector<int> faces;
    const int iface1 = newCells[k][0][0];
    const int iface2 = newCells[k][1][0];
    if (mNodePositions[iface1].x() < cellX[numOldZones + k]) {
      faces.push_back(~iface1);
      faces.push_back(iface2);
    } else {
      faces.push_back(~iface2);
      faces.push_back(iface1);
    }
    unsigned izone = mZones.size();
    mZones.push_back(Zone(*this, izone, faces));
  }

  // Post-conditions.
  ENSURE2(this->valid() == "", this->valid());
}

//------------------------------------------------------------------------------
// Mesh::reconstruct(generators, boundary)
//------------------------------------------------------------------------------
template<>
void
Mesh<Dim<1> >::
reconstructInternal(const vector<Mesh<Dim<1> >::Vector>& localGenerators,
                    const Dim<1>::FacetedVolume& boundary) {
  return this->reconstructInternal(localGenerators, boundary.xmin(), boundary.xmax());
}

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<> const unsigned Mesh<Dim<1> >::minFacesPerZone = 2;
template<> const unsigned Mesh<Dim<1> >::minEdgesPerZone = 2;
template<> const unsigned Mesh<Dim<1> >::minNodesPerZone = 2;
template<> const unsigned Mesh<Dim<1> >::minEdgesPerFace = 1;
template<> const unsigned Mesh<Dim<1> >::minNodesPerFace = 1;

}

//------------------------------------------------------------------------------
// Instantiate the generic mesh non-inlined methods.
//------------------------------------------------------------------------------
#include "Mesh.cc"
template class Spheral::Mesh<Spheral::Dim<1> >;
