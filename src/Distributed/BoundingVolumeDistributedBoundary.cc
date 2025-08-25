//---------------------------------Spheral++----------------------------------//
// BoundingVolumeDistributedBoundary
//
// Build a distributed boundary based on testing for intersecting bounding 
// volumes of domains.
//
// Created by JMO, Tue Jan 19 09:22:37 PST 2010
//----------------------------------------------------------------------------//
#include <mpi.h>

#include "DistributedBoundary.hh"
#include "BoundingVolumeDistributedBoundary.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/NestedGridNeighbor.hh"
#include "Neighbor/GridCellIndex.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Utilities/removeElements.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/nodeBoundingBoxes.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/DBC.hh"
#include "waitAllWithDeadlockDetection.hh"
#include "Communicator.hh"

#include <algorithm>
using std::vector;
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
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
BoundingVolumeDistributedBoundary<Dimension>::
BoundingVolumeDistributedBoundary():
  DistributedBoundary<Dimension>() {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
BoundingVolumeDistributedBoundary<Dimension>::
~BoundingVolumeDistributedBoundary() {
}

//------------------------------------------------------------------------------
// Set the ghost nodes for the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
BoundingVolumeDistributedBoundary<Dimension>::
setAllGhostNodes(DataBase<Dimension>& dataBase) {

  // Clear out the existing communication map for the given database.
  this->reset(dataBase);

  // Start out by determining the bounding volumes for all domains, and then intersecting
  // them to open communication maps of who talks to who.
  buildSendNodes(dataBase);

  // Tell everyone else the send nodes we have for them.
  this->buildReceiveAndGhostNodes(dataBase);

  // Exchange the minimal info we expect for the NodeLists: mass, position, H
  for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    this->updateGhostNodes(**nodeListItr);
  }

  // And that's all.  At this point each domain knows who it it sending nodes to,
  // what nodes to send them, who it is receiving nodes from, and what nodes it
  // be receiving.
}

//------------------------------------------------------------------------------
// Each process should have the set of occupied grid cells for all other domains
// at this point.  Now we use this information to determine which processes this
// domain should be sending nodes to.
//------------------------------------------------------------------------------
template<typename Dimension>
void
BoundingVolumeDistributedBoundary<Dimension>::
buildSendNodes(const DataBase<Dimension>& dataBase) {

  // This processor's ID.
  int procID = this->domainID();
  int numProcs = this->numDomains();
  CHECK(procID < numProcs);

  // Compute the local bounding volumes.
  typedef typename Dimension::ConvexHull ConvexHull;
  vector<ConvexHull> domainNodeBoundingVolume(numProcs), domainSampleBoundingVolume(numProcs);
  globalBoundingVolumes(dataBase, domainNodeBoundingVolume[procID], domainSampleBoundingVolume[procID]);

  // Globally exchange the bounding polyhedra.
  {
    vector<char> localBuffer;
    packElement(domainNodeBoundingVolume[procID], localBuffer);
    packElement(domainSampleBoundingVolume[procID], localBuffer);
    for (int sendProc = 0; sendProc != numProcs; ++sendProc) {
      unsigned bufSize = localBuffer.size();
      MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, Communicator::communicator());
      if (bufSize > 0) {
        vector<char> buffer = localBuffer;
        buffer.resize(bufSize);
        MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, Communicator::communicator());
        vector<char>::const_iterator itr = buffer.begin();
        unpackElement(domainNodeBoundingVolume[sendProc], itr, buffer.end());
        unpackElement(domainSampleBoundingVolume[sendProc], itr, buffer.end());
        CHECK(itr == buffer.end());
      }
    }
  }

  // Compute our node bounding boxes.
  typedef std::pair<Vector, Vector> Box;
  const FieldList<Dimension, Box> nodeSampleBoxes = nodeBoundingBoxes(dataBase);

  // Iterate over all the other domains and check who has bounding volumes that
  // intersect with our own.
  const FieldList<Dimension, Vector> positions = dataBase.globalPosition();
  for (int neighborProc = 0; neighborProc != numProcs; ++neighborProc) {
    if (neighborProc != procID) {
      if (domainSampleBoundingVolume[procID].intersect(domainNodeBoundingVolume[neighborProc]) or   // I see you
          domainSampleBoundingVolume[neighborProc].intersect(domainNodeBoundingVolume[procID])) {   // You see me

        // This domain overlaps ours, so look for any of our nodes who's boxes
        // intersect the other domain.
        int nodeListi = 0;
        for (typename DataBase<Dimension>::ConstNodeListIterator itr = dataBase.nodeListBegin();
             itr != dataBase.nodeListEnd();
             ++itr, ++nodeListi) {
          const int numNodes = (**itr).numNodes();
          vector<int> sendNodes;
          for (int i = 0; i != numNodes; ++i) {
            if (domainNodeBoundingVolume[neighborProc].intersect(nodeSampleBoxes(nodeListi, i)) or  // I see you
                domainSampleBoundingVolume[neighborProc].contains(positions(nodeListi, i))) {       // You see me
              sendNodes.push_back(i);
            }
          }
          if (sendNodes.size() > 0) {
            DomainBoundaryNodes& domainNodes = this->openDomainBoundaryNodes(&(**itr), neighborProc);
            copy(sendNodes.begin(), sendNodes.end(), back_inserter(domainNodes.sendNodes));
          }
        }
      }
    }
  }

  // Remove the duplicates in our send node lists.
  typedef typename DistributedBoundary<Dimension>::DomainBoundaryNodeMap DomainBoundaryNodeMap;
  typedef typename DistributedBoundary<Dimension>::NodeListDomainBoundaryNodeMap NodeListDomainBoundaryNodeMap;
  NodeListDomainBoundaryNodeMap& nodeListDomainBoundaryNodeMap = this->accessNodeListDomainBoundaryNodeMap();
  for (typename NodeListDomainBoundaryNodeMap::iterator itr0 = nodeListDomainBoundaryNodeMap.begin();
       itr0 != nodeListDomainBoundaryNodeMap.end();
       ++itr0) {
    DomainBoundaryNodeMap& domBoundaryNodeMap = itr0->second;
    for (typename DomainBoundaryNodeMap::iterator itr1 = domBoundaryNodeMap.begin();
         itr1 != domBoundaryNodeMap.end();
         ++itr1) {
      vector<size_t>& sendNodes = itr1->second.sendNodes;
      sort(sendNodes.begin(), sendNodes.end());
      sendNodes.erase(unique(sendNodes.begin(), sendNodes.end()), sendNodes.end());
    }
  }

}

//------------------------------------------------------------------------------
// Pack up the positions and H's from a NodeList for communciation.
//------------------------------------------------------------------------------
template<typename Dimension>
void
BoundingVolumeDistributedBoundary<Dimension>::
packNodeListBuffers(const DataBase<Dimension>& dataBase,
                    vector<int>& numNodesPerNodeList,
                    vector<vector<char>>& positionBuffers,
                    vector<vector<char>>& Hbuffers) const {

  for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    numNodesPerNodeList.push_back((**nodeListItr).numInternalNodes());
    positionBuffers.push_back((**nodeListItr).positions().serialize());
    Hbuffers.push_back((**nodeListItr).Hfield().serialize());
  }

  ENSURE(numNodesPerNodeList.size() == dataBase.numNodeLists());
  ENSURE(positionBuffers.size() == numNodesPerNodeList.size());
  ENSURE(Hbuffers.size() == numNodesPerNodeList.size());
}

}

