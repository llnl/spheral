//---------------------------------Spheral++----------------------------------//
// TreeDistributedBoundary -- Implementation of the Distributed Boundary
// condition for use with TreeNeighbor based NodeLists.
//
// Created by JMO, Mon Aug 27 21:57:51 PDT 2001
//----------------------------------------------------------------------------//
#include <mpi.h>

#include "DistributedBoundary.hh"
#include "TreeDistributedBoundary.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/TreeNeighbor.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/removeElements.hh"
#include "Utilities/DBC.hh"
#include "waitAllWithDeadlockDetection.hh"
#include "Communicator.hh"

#include <list>
#include <algorithm>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Singleton instance method.
//------------------------------------------------------------------------------
template<typename Dimension>
TreeDistributedBoundary<Dimension>&
TreeDistributedBoundary<Dimension>::
instance() {
  static TreeDistributedBoundary<Dimension> theInstance;
  return theInstance;
}

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TreeDistributedBoundary<Dimension>::TreeDistributedBoundary():
  DistributedBoundary<Dimension>() {
}

//------------------------------------------------------------------------------
// Set the ghost nodes for the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeDistributedBoundary<Dimension>::
setAllGhostNodes(DataBase<Dimension>& dataBase) {

  // This processor's ID.
  const int procID = this->domainID();
  const int numProcs = this->numDomains();
  CONTRACT_VAR(numProcs);
  CONTRACT_VAR(procID);
  CHECK(procID < numProcs);

  // Clear out the existing communication map for the given database.
  this->reset(dataBase);

  // Get the set of occupied tree cells (by key)
  const auto* neighborPtr = this->getTreeNeighborPtr(dataBase);
  auto localTree = neighborPtr->occupiedCells();

  // Use per domain trees to build the local send nodes needed for each domain.
  this->buildSendNodes(dataBase, localTree);

  // Tell everyone else the send nodes we have for them.
  this->buildReceiveAndGhostNodes(dataBase);

  // Exchange the minimal info we expect for the NodeLists: mass, position, H
  for (auto nodeListItr = dataBase.nodeListBegin(); nodeListItr != dataBase.nodeListEnd(); ++nodeListItr) {
    this->updateGhostNodes(**nodeListItr);
  }

  // And that's all.  At this point each domain knows who it it sending nodes to,
  // what nodes to send them, who it is receiving nodes from, and what nodes it will
  // be receiving.
}

//------------------------------------------------------------------------------
// Get the TreeNeighbor associated with the given NodeList, or throw up if
// it's not there.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TreeNeighbor<Dimension>*
TreeDistributedBoundary<Dimension>::
getTreeNeighborPtr(const DataBase<Dimension>& dataBase) const {
  const auto* result = std::get_if<TreeNeighbor<Dimension>>(&dataBase.neighbor());
  VERIFY2(result != nullptr, "TreeDistributedBoundary ERROR : unable to extract TreeNeighbor from DataBase");
  return result;
}

//------------------------------------------------------------------------------
// Each process should have the set of occupied grid cells for all other domains
// at this point.  Now we use this information to determine which processes this
// domain should be sending nodes to.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeDistributedBoundary<Dimension>::
buildSendNodes(const DataBase<Dimension>& dataBase,
               const vector<vector<typename TreeDistributedBoundary<Dimension>::CellKey>>& localTree) {

  // This processor's ID.
  int procID = this->domainID();
  int numProcs = this->numDomains();
  CHECK(procID < numProcs);

  // Get the state from the DataBase
  const auto* neighborPtr = this->getTreeNeighborPtr(dataBase);
  const auto& nodeListPtrs = dataBase.nodeListPtrs();
  const auto& positions = dataBase.globalPosition();
  const auto& extents = neighborPtr->nodeExtentView();

  // Compute the local H inverse for our box culling of send nodes.
  FieldList<Dimension, SymTensor> Hinverse = dataBase.newGlobalFieldList(SymTensor::zero(), "H inverse");
  dataBase.globalHinverse(Hinverse);

  // Serialize our local tree for broadcast.
  vector<char> localBuffer;
  packElement(localTree, localBuffer);

  // We also pack the min/max node extents.
  // This is used to cull the nodes we're going to send to other domains.
  double radiusNodesLocal, radiusSampleLocal;
  Vector centroidLocal, xminNodesLocal, xmaxNodesLocal, xminSampleLocal, xmaxSampleLocal;
  dataBase.localSamplingBoundingVolume(centroidLocal, 
                                       radiusNodesLocal,
                                       radiusSampleLocal,
                                       xminNodesLocal,
                                       xmaxNodesLocal,
                                       xminSampleLocal,
                                       xmaxSampleLocal);
  packElement(centroidLocal, localBuffer);
  packElement(radiusNodesLocal, localBuffer);
  packElement(radiusSampleLocal, localBuffer);
  packElement(xminNodesLocal, localBuffer);
  packElement(xmaxNodesLocal, localBuffer);
  packElement(xminSampleLocal, localBuffer);
  packElement(xmaxSampleLocal, localBuffer);
  unsigned localBufSize = localBuffer.size();
  
  // Fire off a storm of non-blocking messages to all other domains with our packed occupancy info.
  vector<MPI_Request> sendRequests;
  for (auto sendProc = 0; sendProc < numProcs; ++sendProc) {
    if (sendProc != procID) {
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&localBufSize, 1, MPI_UNSIGNED, sendProc, procID, Communicator::communicator(), &sendRequests.back());
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&localBuffer[0], localBufSize, MPI_CHAR, sendProc, numProcs + procID, Communicator::communicator(), &sendRequests.back());
    }
  }
  CHECK((int)sendRequests.size() == 2*(numProcs -1 ));

  // Next walk the other domains, get their info, and figure out what nodes we need to send them.
  for (auto recvProc = 0; recvProc < numProcs; ++recvProc) {
    if (recvProc != procID) {

      // Get the neighbor serialized info.
      MPI_Status status1, status2;
      unsigned bufSize;
      MPI_Recv(&bufSize, 1, MPI_UNSIGNED, recvProc, recvProc, Communicator::communicator(), &status1);
      vector<char> buffer(bufSize);
      MPI_Recv(&buffer[0], bufSize, MPI_CHAR, recvProc, numProcs + recvProc, Communicator::communicator(), &status2);

      // Unpack the neighbor tree and bounds.
      vector<vector<CellKey>> otherTree;
      double radiusNodesOther, radiusSampleOther;
      Vector centroidOther, xminNodesOther, xmaxNodesOther, xminSampleOther, xmaxSampleOther;
      vector<char>::const_iterator bufItr = buffer.begin();
      unpackElement(otherTree, bufItr, buffer.end());
      unpackElement(centroidOther, bufItr, buffer.end());
      unpackElement(radiusNodesOther, bufItr, buffer.end());
      unpackElement(radiusSampleOther, bufItr, buffer.end());
      unpackElement(xminNodesOther, bufItr, buffer.end());
      unpackElement(xmaxNodesOther, bufItr, buffer.end());
      unpackElement(xminSampleOther, bufItr, buffer.end());
      unpackElement(xmaxSampleOther, bufItr, buffer.end());
      CHECK(bufItr == buffer.end());

      // Walk the local tree and see if we interact with the remote tree
      for (auto klevel = 0u; klevel < otherTree.size(); ++klevel) {
        for (auto kcell = 0u; kcell < otherTree[klevel].size(); ++kcell) {

          // Find the master/coarse neighbor info for this tree level/cell.
          std::vector<NodeID> masterList, coarseNeighbors;
          neighborPtr->setTreeMasterList(klevel, otherTree[klevel][kcell], masterList, coarseNeighbors, false);

          // Copy the coarse neighbor set as our send nodes.
          for (const auto& nodeID: coarseNeighbors) {
            CHECK(nodeID.first < nodeListPtrs.size());
            auto& domainNodes = this->openDomainBoundaryNodes(nodeListPtrs[nodeID.first], recvProc);
            domainNodes.sendNodes.push_back(nodeID.second);
          }
        }
      }

      // Check our sendNodes for this process
      for (auto [k, nodeListPtr]: enumerate(nodeListPtrs)) {
        // Did we identify any send nodes for this NodeList?
        if (this->nodeListSharedWithDomain(*nodeListPtr, recvProc)) {

          // Remove any duplicate nodes that were created.
          auto& sendNodes = this->accessDomainBoundaryNodes(*nodeListPtr, recvProc).sendNodes;
	  CHECK(sendNodes.size() > 0);
          sort(sendNodes.begin(), sendNodes.end());
          sendNodes.erase(std::unique(sendNodes.begin(), sendNodes.end()), sendNodes.end());
	  CHECK(sendNodes.size() > 0);

          // Reject nodes based on the "box culling" domain of influence logic.
          std::vector<int> indicesToKill;
          for (auto [kk, i]: enumerate(sendNodes)) {
            const auto& xi = positions(k,i);
            const auto& extenti = extents(k,i);
            const auto  dr = xi - centroidOther;
            const auto  drMag = dr.magnitude();
            const auto  drUnit = dr.unitVector();
            const auto  xmini = xi - extenti;
            const auto  xmaxi = xi + extenti;
            const auto  hi = (Hinverse(k,i)*drUnit).magnitude();
            const bool  interacting = (drMag - 2.0*hi < radiusNodesOther or  // I see you
                                       drMag < radiusSampleOther or          // You see me
                                       testBoxIntersection(xmini, xmaxi, xminNodesOther, xmaxNodesOther) or // I see you
                                       testBoxIntersection(xi, xi, xminSampleOther, xmaxSampleOther));      // You see me
            if (not interacting) indicesToKill.push_back(kk);
          }

          // Are we killing all of the send nodes?
          if (indicesToKill.size() == sendNodes.size()) {
            this->removeDomainBoundaryNodes(nodeListPtr, recvProc);
          } else {
            removeElements(sendNodes, indicesToKill);
          }
          // cerr << "  --> Sending to " << recvProc << " : ";
          // std::copy(sendNodes.begin(), sendNodes.end(), std::ostream_iterator<int>(std::cerr, " "));
          // cerr << endl;
        }
      }
    }
  }

  // Wait until all our send buffers have been received.
  if (not sendRequests.empty()) {
    vector<MPI_Status> status(sendRequests.size());
    MPI_Waitall(sendRequests.size(), &sendRequests[0], &status[0]);
  }
}

}

