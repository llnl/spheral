//---------------------------------Spheral++----------------------------------//
// RedistributeNodes -- An abstract base class for methods that repartition
// the Spheral++ NodeLists among domains.
//
// Created by JMO, Tue Feb  4 14:23:11 PST 2003
//----------------------------------------------------------------------------//
#include "RedistributeNodes.hh"
#include "Utilities/DomainNode.hh"
#include "DataBase/DataBase.hh"
#include "NodeList/NodeList.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/packElement.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/RedistributionRegistrar.hh"
#include "Communicator.hh"

#include "Utilities/DBC.hh"

#include <mpi.h>

#include <algorithm>
#include <vector>
#include <list>
#include <sstream>
using std::vector;
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

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
template<typename Dimension>
RedistributeNodes<Dimension>::
RedistributeNodes():
  mComputeWork(false) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
RedistributeNodes<Dimension>::
~RedistributeNodes() {
}

//------------------------------------------------------------------------------
// Get the total number of nodes in the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
int
RedistributeNodes<Dimension>::
numGlobalNodes(const DataBase<Dimension>& dataBase) const {

  typedef typename DataBase<Dimension>::ConstNodeListIterator ItrType;
  int result = 0;
  for (ItrType itr = dataBase.nodeListBegin();
       itr < dataBase.nodeListEnd();
       ++itr) {
    result += Spheral::numGlobalNodes(**itr);
  }
  return result;
}

//------------------------------------------------------------------------------
// Calculate the current domain decomposition, and return it as a set of 
// DomainNode identifiers.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<DomainNode<Dimension>>
RedistributeNodes<Dimension>::
currentDomainDecomposition(const DataBase<Dimension>& dataBase,
                           const FieldList<Dimension, size_t>& globalNodeIDs) const {
  FieldList<Dimension, Scalar> dummyWork = dataBase.newGlobalFieldList(Scalar());
  return currentDomainDecomposition(dataBase, globalNodeIDs, dummyWork);
}

//------------------------------------------------------------------------------
// Calculate the current domain decomposition, and return it as a set of 
// DomainNode identifiers.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<DomainNode<Dimension> >
RedistributeNodes<Dimension>::
currentDomainDecomposition(const DataBase<Dimension>& dataBase,
                           const FieldList<Dimension, size_t>& globalNodeIDs,
                           const FieldList<Dimension, Scalar>& workPerNode) const {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(dataBase.numNodeLists() == globalNodeIDs.numFields());
  REQUIRE(dataBase.numNodeLists() == workPerNode.numFields());
  for (auto* fptr: globalNodeIDs) {
    CONTRACT_VAR(fptr) ;
    REQUIRE(dataBase.haveNodeList(fptr->nodeList()));
  }
  for (auto* fptr: workPerNode) {
    CONTRACT_VAR(fptr) ;
    REQUIRE(dataBase.haveNodeList(fptr->nodeList()));
  }
  END_CONTRACT_SCOPE

  // Prepare the result.
  vector<DomainNode<Dimension>> result;

  // Loop over each NodeList in the DataBase.
  int nodeListID = 0;
  const int proc = domainID();
  int offset = 0;
  for (auto* nodeListPtr: dataBase.nodeListPtrs()) {
    const auto& nodeList = *nodeListPtr;
    const auto& globalNodeListIDs = **(globalNodeIDs.fieldForNodeList(nodeList));
    const auto& work = **(workPerNode.fieldForNodeList(nodeList));

    // Loop over the nodes in the NodeList, and add their info to the 
    // result.
    result.reserve(result.size() + nodeList.numInternalNodes());
    for (auto localID = 0u; localID != nodeList.numInternalNodes(); ++localID) {
      result.push_back(DomainNode<Dimension>());
      result.back().localNodeID = localID;
      result.back().uniqueLocalNodeID = localID + offset;
      result.back().globalNodeID = globalNodeListIDs(localID);
      result.back().nodeListID = nodeListID;
      result.back().domainID = proc;
      result.back().work = work(localID);
      result.back().position = nodeList.positions()(localID);
    }
    offset += nodeList.numInternalNodes();
    ++nodeListID;
  }
  CHECK(nodeListID == (int)dataBase.numNodeLists());

  ENSURE(validDomainDecomposition(result, dataBase));
  return result;
}

//------------------------------------------------------------------------------
// Given a desired domain decomposition, reassign the nodes across domains and
// create it!
//------------------------------------------------------------------------------
template<typename Dimension>
void
RedistributeNodes<Dimension>::
enforceDomainDecomposition(const vector<DomainNode<Dimension> >& nodeDistribution,
                           DataBase<Dimension>& dataBase) const {

  // Notify everyone before we start redistributing.
  RedistributionRegistrar::instance().preRedistributionNotifications();

  REQUIRE(validDomainDecomposition(nodeDistribution, dataBase));
  typedef typename DataBase<Dimension>::ConstNodeListIterator NodeListIterator;

  const int numNodeLists = dataBase.numNodeLists();
  const int proc = domainID();
  const int numProcs = numDomains();

  // Post receives for how many nodes we're getting from each domain.
  vector< vector<int> > numRecvNodes;
  numRecvNodes.reserve(numProcs);
  for (int i = 0; i != numProcs; ++i) numRecvNodes.push_back(vector<int>(numNodeLists, 0));
  CHECK((int)numRecvNodes.size() == numProcs);
  vector<MPI_Request> numRecvNodeRequests;
  numRecvNodeRequests.reserve(numProcs - 1);
  for (int recvProc = 0; recvProc != numProcs; ++recvProc) {
    CHECK(recvProc < (int)numRecvNodes.size());
    if (recvProc != proc) {
      numRecvNodeRequests.push_back(MPI_Request());
      MPI_Irecv(&(*numRecvNodes[recvProc].begin()), numNodeLists, MPI_INT, recvProc, 0, Communicator::communicator(), &(numRecvNodeRequests.back()));
    }
  }
  CHECK((int)numRecvNodeRequests.size() == numProcs - 1);

  // Find the nodes on this domain that have to be reassigned to new domains.
  vector< vector< vector<size_t> > > sendNodes(numProcs); // [sendDomain][nodeList][node]
  for (int i = 0; i != numProcs; ++i) sendNodes[i].resize(numNodeLists);
  for (typename vector<DomainNode<Dimension> >::const_iterator distItr = nodeDistribution.begin();
       distItr != nodeDistribution.end();
       ++distItr) {
    if (distItr->domainID != proc) {
      CHECK(distItr->domainID >= 0 && distItr->domainID < numProcs);
      CHECK(distItr->nodeListID < numNodeLists);
      CHECK(distItr->domainID < (int)sendNodes.size());
      CHECK(distItr->nodeListID < (int)sendNodes[distItr->domainID].size());
      sendNodes[distItr->domainID][distItr->nodeListID].push_back(distItr->localNodeID);
    }
  }

  // Tell everyone the number of nodes we're sending them.
  int numSendDomains = 0;
  list< vector<int> > numSendNodes;
  vector<MPI_Request> numSendNodeRequests;
  numSendNodeRequests.reserve(numProcs - 1);
  vector<int> totalNumSendNodes(numProcs, 0);
  for (int sendProc = 0; sendProc != numProcs; ++sendProc) {
    if (sendProc != proc) {
      numSendNodes.push_back(vector<int>());
      vector<int>& x = numSendNodes.back();
      x.reserve(numNodeLists);
      for (int nodeListID = 0; nodeListID != numNodeLists; ++nodeListID) {
        x.push_back(sendNodes[sendProc][nodeListID].size());
        totalNumSendNodes[sendProc] += x.back();
      }
      CHECK((int)x.size() == numNodeLists);
      numSendNodeRequests.push_back(MPI_Request());
      MPI_Isend(&(*x.begin()), numNodeLists, MPI_INT, sendProc, 0, Communicator::communicator(), &(numSendNodeRequests.back()));
      if (totalNumSendNodes[sendProc] > 0) ++numSendDomains;
    }
  }
  CHECK((int)numSendNodes.size() == numProcs - 1)
  CHECK((int)numSendNodeRequests.size() == numProcs - 1);
  CHECK(numSendDomains <= numProcs - 1);
  CHECK(totalNumSendNodes[proc] == 0);

  // Pack up the field info for the nodes we're sending.
  vector< vector< list< vector<char> > > > sendBuffers(numProcs);
  for (int sendProc = 0; sendProc != numProcs; ++sendProc) {
    CHECK(sendProc < (int)sendBuffers.size());
    sendBuffers[sendProc].resize(numNodeLists);
    int nodeListID = 0;
    for (NodeListIterator nodeListItr = dataBase.nodeListBegin();
         nodeListItr != dataBase.nodeListEnd();
         ++nodeListItr, ++nodeListID) {
      const NodeList<Dimension>& nodeList = **nodeListItr;
      CHECK(sendProc < (int)sendNodes.size());
      CHECK(nodeListID < (int)sendNodes[sendProc].size());
      const int numSendNodes = sendNodes[sendProc][nodeListID].size();
      if (numSendNodes > 0) {
        CHECK(sendProc < (int)sendBuffers.size());
        CHECK(nodeListID < (int)sendBuffers[sendProc].size());
        sendBuffers[sendProc][nodeListID] = nodeList.packNodeFieldValues(sendNodes[sendProc][nodeListID]);
      }
    }
    CHECK(nodeListID == numNodeLists);
  }

  // Wait until we know who is sending us nodes.
  if (not numRecvNodeRequests.empty()) {
    vector<MPI_Status> recvStatus(numRecvNodeRequests.size());
    MPI_Waitall(numRecvNodeRequests.size(), &(*numRecvNodeRequests.begin()), &(*recvStatus.begin()));
  }

  // Sum the total nodes we're receiving from each domain.
  int numRecvDomains = 0;
  vector<int> totalNumRecvNodes(numProcs, 0);
  for (int recvProc = 0; recvProc != numProcs; ++recvProc) {
    if (recvProc != proc) {
      CHECK(recvProc < (int)totalNumRecvNodes.size());
      CHECK(recvProc < (int)numRecvNodes.size());
      CHECK((int)numRecvNodes[recvProc].size() == numNodeLists);
      for (vector<int>::const_iterator itr = numRecvNodes[recvProc].begin();
           itr != numRecvNodes[recvProc].end();
           ++itr) totalNumRecvNodes[recvProc] += *itr;
      if (totalNumRecvNodes[recvProc] > 0) ++numRecvDomains;
    }
  }
  CHECK(numRecvDomains >= 0);
  CHECK(totalNumRecvNodes[proc] == 0);

  // Post receives for the sizes of the encoded field buffers we'll be receiving.
  int numRecvBuffers = 0;
  list< list< vector<int> > > recvBufferSizes;
  vector<MPI_Request> recvBufSizeRequests;
  recvBufSizeRequests.reserve(numRecvDomains*numNodeLists);
  for (int recvProc = 0; recvProc != numProcs; ++recvProc) {
    if (totalNumRecvNodes[recvProc] > 0) {
      recvBufferSizes.push_back(list< vector<int> >());
      list< vector<int> >& x = recvBufferSizes.back();
      int nodeListID = 0;
      for (NodeListIterator itr = dataBase.nodeListBegin();
           itr != dataBase.nodeListEnd();
           ++itr, ++nodeListID) {
        if (numRecvNodes[recvProc][nodeListID] > 0) {
          const int numFields = (*itr)->numFields();
          x.push_back(vector<int>(numFields));
          recvBufSizeRequests.push_back(MPI_Request());
          MPI_Irecv(&(*x.back().begin()), numFields, MPI_INT, recvProc, 1 + nodeListID, Communicator::communicator(), &(recvBufSizeRequests.back()));
          numRecvBuffers += numFields;
        }
      }
      CHECK((int)x.size() <= numNodeLists);
      CHECK(nodeListID == numNodeLists);
    }
  }
  CHECK((int)recvBufferSizes.size() == numRecvDomains);
  CHECK((int)recvBufSizeRequests.size() <= numRecvDomains*numNodeLists);
  CHECK(numRecvBuffers >= numRecvDomains);

  // Send the sizes of the our encoded field buffers.
  int numSendBuffers = 0;
  list< list< vector<int> > > sendBufferSizes;
  vector<MPI_Request> sendBufSizeRequests;
  sendBufSizeRequests.reserve(numSendDomains*numNodeLists);
  for (int sendProc = 0; sendProc != numProcs; ++sendProc) {
    if (totalNumSendNodes[sendProc] > 0) {
      sendBufferSizes.push_back(list< vector<int> >());
      list< vector<int> >& x = sendBufferSizes.back();
      for (int nodeListID = 0; nodeListID != numNodeLists; ++nodeListID) {
        if (sendNodes[sendProc][nodeListID].size() > 0) {
          const list< vector<char> >& sendBufs = sendBuffers[sendProc][nodeListID];
          x.push_back(vector<int>());
          for (list< vector<char> >::const_iterator itr = sendBufs.begin();
               itr != sendBufs.end();
               ++itr) {
            x.back().push_back(itr->size());
          }
          numSendBuffers += sendBufs.size();
          CHECK(x.back().size() == sendBufs.size());
          sendBufSizeRequests.push_back(MPI_Request());
          MPI_Isend(&(*x.back().begin()), x.back().size(), MPI_INT, sendProc, 1 + nodeListID, Communicator::communicator(), &(sendBufSizeRequests.back()));
        }
      }
      CHECK((int)x.size() <= numNodeLists);
    }
  }
  CHECK((int)sendBufferSizes.size() == numSendDomains);
  CHECK((int)sendBufSizeRequests.size() <= numSendDomains*numNodeLists);
  CHECK(numSendBuffers >= numSendDomains);
      
  // Determine the maximum number of fields defined on a NodeList.
  size_t maxNumFields = 0u;
  for (auto* nodeListPtr: dataBase.nodeListPtrs()) {
    maxNumFields = std::max(maxNumFields, nodeListPtr->numFields());

    // This is a somewhat expensive contract because of the all reduce,
    // but it is critical all processors agree about the fields defined
    // on each NodeList.
    BEGIN_CONTRACT_SCOPE
    {
      int localNumFields = nodeListPtr->numFields();
      int globalNumFields;
      MPI_Allreduce(&localNumFields, &globalNumFields, 1, MPI_INT, MPI_MAX, Communicator::communicator());
      CHECK(localNumFields == globalNumFields);
    }
    END_CONTRACT_SCOPE
  }

  // Wait until we know the sizes of the encoded field buffers we're receiving.
  if (not recvBufSizeRequests.empty()) {
    vector<MPI_Status> recvStatus(recvBufSizeRequests.size());
    MPI_Waitall(recvBufSizeRequests.size(), &(*recvBufSizeRequests.begin()), &(*recvStatus.begin()));
  }

  // Prepare to receive the encoded field buffers.
  list< list< list< vector<char> > > > fieldBuffers;
  list< list< vector<int> > >::const_iterator outerBufSizeItr = recvBufferSizes.begin();
  vector<MPI_Request> recvBufferRequests;
  recvBufferRequests.reserve(numRecvBuffers);
  for (int recvProc = 0; recvProc != numProcs; ++recvProc) {
    if (totalNumRecvNodes[recvProc] > 0) {
      fieldBuffers.push_back(list< list< vector<char> > >());
      CHECK(outerBufSizeItr != recvBufferSizes.end());
      CHECK((int)outerBufSizeItr->size() <= numNodeLists);
      list< vector<int> >::const_iterator bufSizeItr = outerBufSizeItr->begin();
      for (int nodeListID = 0; nodeListID != numNodeLists; ++nodeListID) {
        if (numRecvNodes[recvProc][nodeListID] > 0) {
          CHECK(bufSizeItr != outerBufSizeItr->end());
          fieldBuffers.back().push_back(list< vector<char> >());
          list< vector<char> >& bufs = fieldBuffers.back().back();
          const vector<int>& bufSizes = *bufSizeItr;
          size_t i = 0;
          for (vector<int>::const_iterator itr = bufSizes.begin();
               itr != bufSizes.end();
               ++itr, ++i) {
            bufs.push_back(vector<char>(*itr));
            recvBufferRequests.push_back(MPI_Request());
            MPI_Irecv(&(*bufs.back().begin()), *itr, MPI_CHAR, recvProc, (2 + numNodeLists) + maxNumFields*nodeListID + i, Communicator::communicator(), &(recvBufferRequests.back()));
          }
          ++bufSizeItr;
          CHECK(i <= maxNumFields);
        }
      }
      CHECK(bufSizeItr == outerBufSizeItr->end());
      ++outerBufSizeItr;
    }
  }
  CHECK((int)fieldBuffers.size() == numRecvDomains);
  CHECK(outerBufSizeItr == recvBufferSizes.end());
  CHECK((int)recvBufferRequests.size() == numRecvBuffers);

  // Send our encoded buffers.
  vector<MPI_Request> sendBufferRequests;
  sendBufferRequests.reserve(numSendBuffers);
  for (int sendProc = 0; sendProc != numProcs; ++sendProc) {
    if (totalNumSendNodes[sendProc] > 0) {
      CHECK(sendProc < (int)sendBuffers.size());
      CHECK((int)sendBuffers[sendProc].size() == numNodeLists);
      for (int nodeListID = 0; nodeListID != numNodeLists; ++nodeListID) {
        if (sendNodes[sendProc][nodeListID].size() > 0) {
          list< vector<char> >& bufs = sendBuffers[sendProc][nodeListID];
          size_t i = 0;
          for (list< vector<char> >::iterator itr = bufs.begin();
               itr != bufs.end();
               ++itr, ++i) {
            vector<char>& buf = *itr;
            sendBufferRequests.push_back(MPI_Request());
            MPI_Isend(&(*buf.begin()), buf.size(), MPI_CHAR, sendProc, (2 + numNodeLists) + maxNumFields*nodeListID + i, Communicator::communicator(), &(sendBufferRequests.back()));
          }
          CHECK(i <= maxNumFields);
        }
      }
    }
  }
  CHECK((int)sendBufferRequests.size() == numSendBuffers);

  // Wait until we've received the packed field buffers.
  if (not recvBufferRequests.empty()) {
    vector<MPI_Status> recvStatus(recvBufferRequests.size());
    MPI_Waitall(recvBufferRequests.size(), &(*recvBufferRequests.begin()), &(*recvStatus.begin()));
  }

  // OK, now we're done with the communication part of this scheme.  Next we need to 
  // delete the nodes (and their associated Field elements) that we've transferred
  // to other domains.
  {
    int nodeListID = 0;
    for (NodeListIterator nodeListItr = dataBase.nodeListBegin();
         nodeListItr != dataBase.nodeListEnd();
         ++nodeListItr, ++nodeListID) {
      NodeList<Dimension>& nodeList = **nodeListItr;

      // Build the complete list of local IDs we're removing from this NodeList.
      // You have to do this, 'cause if you delete just some nodes from the 
      // NodeList, then the remaining sendNode localIDs are not valid.
      vector<size_t> deleteNodes;
      for (int sendProc = 0; sendProc != numProcs; ++sendProc) {
        CHECK(sendProc < (int)sendNodes.size());
        CHECK(nodeListID < (int)sendNodes[sendProc].size());
        deleteNodes.reserve(deleteNodes.size() + 
                            sendNodes[sendProc][nodeListID].size());
        for (auto i: sendNodes[sendProc][nodeListID]) deleteNodes.push_back(i);
      }

      // Delete these nodes from the NodeList.
      nodeList.deleteNodes(deleteNodes);
    }
    CHECK(nodeListID == numNodeLists);
  }

  // Now unpack the new nodes and Field values.
  list< list< list< vector<char> > > >::const_iterator procBufItr = fieldBuffers.begin();
  for (int recvProc = 0; recvProc != numProcs; ++recvProc) {
    CHECK(recvProc < (int)totalNumRecvNodes.size());
    if (totalNumRecvNodes[recvProc] > 0) {
      CHECK(procBufItr != fieldBuffers.end());
      int nodeListID = 0;
      list< list< vector<char> > >::const_iterator nodeListBufItr = procBufItr->begin();
      for (NodeListIterator nodeListItr = dataBase.nodeListBegin();
           nodeListItr != dataBase.nodeListEnd();
           ++nodeListItr, ++nodeListID) {
        const int numNewNodes = numRecvNodes[recvProc][nodeListID];
        if (numNewNodes > 0) {
          CHECK(nodeListBufItr != procBufItr->end());
          const list< vector<char> >& bufs = *nodeListBufItr;
          (*nodeListItr)->appendInternalNodes(numNewNodes, bufs);
          (*nodeListItr)->neighbor().updateNodes();
          ++nodeListBufItr;
        }
      }
      CHECK(nodeListBufItr == procBufItr->end());
      ++procBufItr;
    }
  }
  CHECK(procBufItr == fieldBuffers.end());

  // Wait until all our sends are completed.
  if (not numSendNodeRequests.empty()) {
    vector<MPI_Status> sendStatus(numSendNodeRequests.size());
    MPI_Waitall(numSendNodeRequests.size(), &(*numSendNodeRequests.begin()), &(*sendStatus.begin()));
  }
  if (not sendBufSizeRequests.empty()) {
    vector<MPI_Status> sendStatus(sendBufSizeRequests.size());
    MPI_Waitall(sendBufSizeRequests.size(), &(*sendBufSizeRequests.begin()), &(*sendStatus.begin()));
  }
  if (not sendBufferRequests.empty()) {
    vector<MPI_Status> sendStatus(sendBufferRequests.size());
    MPI_Waitall(sendBufferRequests.size(), &(*sendBufferRequests.begin()), &(*sendStatus.begin()));
  }

  // Notify everyone that the nodes have just been shuffled around.
  RedistributionRegistrar::instance().broadcastRedistributionNotifications();
}

//------------------------------------------------------------------------------
// Test that the given domain decomposition is valid (all nodes accounted 
// for, once and only once, etc.).
//------------------------------------------------------------------------------
template<typename Dimension>
bool
RedistributeNodes<Dimension>::
validDomainDecomposition(const vector<DomainNode<Dimension> >& nodeDistribution,
                         const DataBase<Dimension>& dataBase) const {

  // Parallel domain info.
  const int proc = domainID();
  const int numProcs = numDomains();

  // Data structure for checking the local node IDs.
  const int numNodeLists = dataBase.numNodeLists();
  vector< vector<bool> > localIDtaken(numNodeLists);
  int nodeListID = 0;
  for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = 
         dataBase.nodeListBegin();
       nodeListItr < dataBase.nodeListEnd();
       ++nodeListItr, ++nodeListID) {
    CHECK(nodeListID >= 0 && nodeListID < (int)localIDtaken.size());
    localIDtaken[nodeListID].resize((*nodeListItr)->numInternalNodes(), false);
  }

  // Get the maximum allowed value for global IDs.
  const int maxGlobalNodeID = numGlobalNodes(dataBase);

  // We'll build up a vector of just the global node ID assignments.
  vector<int> globalNodeIDs;
  globalNodeIDs.reserve(nodeDistribution.size());

  // Go over every domain node assignment.
  bool valid = true;
  typename vector<DomainNode<Dimension> >::const_iterator nodeDistItr = nodeDistribution.begin();
  while (nodeDistItr < nodeDistribution.end() && valid) {
    const DomainNode<Dimension>& domainNode = *nodeDistItr;
    const int localNodeID = domainNode.localNodeID;
    const int globalNodeID = domainNode.globalNodeID;
    const int nodeListID = domainNode.nodeListID;
    const int domain = domainNode.domainID;

    // Check that the nodeListID is sensible.
    if (!(nodeListID >= 0 && nodeListID < numNodeLists)) valid = false;
    CHECK(valid);

    // Check that this local ID is not already used.
    CHECK(nodeListID >= 0 && nodeListID < numNodeLists);
    CHECK(localNodeID >= 0 && localNodeID < (int)localIDtaken[nodeListID].size());
    if (localIDtaken[nodeListID][localNodeID] == true) {
      valid = false;
    } else {
      localIDtaken[nodeListID][localNodeID] = true;
    }
    CHECK(valid);

    // Check if this global node ID is already allocated on this domain,
    // and insert it into the list of global node IDs for this domain.
    CHECK(globalNodeID >= 0 && globalNodeID < maxGlobalNodeID);
    valid = valid && globalNodeID >= 0 && globalNodeID < maxGlobalNodeID;
    CHECK(valid);
    if (count(globalNodeIDs.begin(), globalNodeIDs.end(), globalNodeID) > 0)
      valid = false;
    globalNodeIDs.push_back(globalNodeID);
    if (!valid) {
      cerr << globalNodeID << " : ";
      for (vector<int>::const_iterator itr = globalNodeIDs.begin();
           itr < globalNodeIDs.end();
           ++itr) cerr << *itr << " ";
      cerr << endl;
    }
    CHECK(valid);

    // Check that the assigned domain is valid.
    valid = valid && (domain >= 0 && domain < numProcs);
    CHECK(valid);

    ++nodeDistItr;
  }

  // Check that all local IDs are accounted for.
  for (vector< vector<bool> >::const_iterator itr = localIDtaken.begin();
       itr < localIDtaken.end();
       ++itr) {
    valid = valid && (count(itr->begin(), itr->end(), false) == 0);
  }
  CHECK(valid);

  // Have all processors agree if we're valid at this point.  If not, go ahead
  // and exit.
  int iValid = 0;
  if (valid) iValid = 1;
  int globalValid;
  MPI_Allreduce(&iValid, &globalValid, 1, MPI_INT, MPI_MIN, Communicator::communicator());
  if (globalValid == 0) {
    return false;
  }

  // Now for the really expensive check.  Make sure that the global IDs 
  // are unique across all processors.
  for (int checkProc = 0; checkProc < numProcs - 1; ++checkProc) {

    if (proc == checkProc) {
      // If we're the processor being checked, then send our list of global
      // Ids to each of the higher processors.
      int numNodes = globalNodeIDs.size();
      for (int sendProc = checkProc + 1; sendProc < numProcs; ++sendProc) {
        MPI_Send(&numNodes, 1, MPI_INT, sendProc, 5, Communicator::communicator());
        MPI_Send(&(*globalNodeIDs.begin()), numNodes, MPI_INT, sendProc,
                 6, Communicator::communicator());
      }

      // Now wait for the responses from each higher processor.
      for (int recvProc = checkProc + 1; recvProc < numProcs; ++recvProc) {
        MPI_Status status;
        int iValid;
        MPI_Recv(&iValid, 1, MPI_INT, recvProc, 7, Communicator::communicator(),
                 &status);
        valid = valid && iValid == 1;
      }
      CHECK(valid);

    } else if (proc > checkProc) {
      // If we're one of the processors greater than the processor being
      // checked, we'll receive the set of global nodes from the check 
      // processor.
      MPI_Status status;
      int numCheckNodes;
      MPI_Recv(&numCheckNodes, 1, MPI_INT, checkProc, 5,
               Communicator::communicator(), &status);
      CHECK(numCheckNodes >= 0);
      vector<int> checkGlobalNodeIDs(numCheckNodes);
      MPI_Recv(&(*checkGlobalNodeIDs.begin()), numCheckNodes, MPI_INT,
               checkProc, 6, Communicator::communicator(), &status);

      // Check to see if any of the global IDs assigned to the processor
      // we're checking is also assigned to this one.
      vector<int>::const_iterator globalCheckItr = checkGlobalNodeIDs.begin();
      while (globalCheckItr < checkGlobalNodeIDs.end() && valid) {
        valid = find(globalNodeIDs.begin(), globalNodeIDs.end(), 
                     *globalCheckItr) == globalNodeIDs.end();
        ++globalCheckItr;
      }
      CHECK(valid);

      // Send the verdict back to the check processor.
      int iValid = 0;
      if (valid) iValid = 1;
      MPI_Send(&iValid, 1, MPI_INT, checkProc, 7, Communicator::communicator());
      CHECK(valid);
    }

    // Now the check processor broadcasts the result of it's check back 
    // to everyone.
    int iValid = 0;
    if (valid) iValid = 1;
    MPI_Bcast(&iValid, 1, MPI_INT, checkProc, Communicator::communicator());
    if (iValid == 1) {
      valid = true;
    } else {
      valid = false;
    }

    // If we not valid, then finish and return false.
    if (!valid) return false;
  }

  // If we got all the way here, it must be valid.
  CHECK(valid);
  return valid;
}

//------------------------------------------------------------------------------
// Assign unique global integer IDs to each node on every domain in the given
// DataBase.  We assign globalIDs so that they are consecutive on a domain,
// increasing with domain.
// Note this is different than what the globalNodeIDs function produces only
// because that one does unique IDs by NodeList.
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Compute the work per node on the local processor.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
RedistributeNodes<Dimension>::
workPerNode(const DataBase<Dimension>& dataBase,
            const double /*Hextent*/) const {

  // Prepare the result.
  FieldList<Dimension, Scalar> result = dataBase.newGlobalFieldList(Scalar(), "work per node");

  // We have two choices for the work: we can use the estimate computed during 
  // the last integration step, or simply use the number of neighbors.
  if (mComputeWork) {

    // The work per node is just the number of neighbors.
    dataBase.updateConnectivityMap(false, false, false);
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
    const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
    for (auto iNodeList = 0u; iNodeList != nodeLists.size(); ++iNodeList) {
      const NodeList<Dimension>* nodeListPtr = nodeLists[iNodeList];
      for (auto i = 0u; i != nodeListPtr->numInternalNodes(); ++i) {
        result(iNodeList, i) = connectivityMap.numNeighborsForNode(nodeListPtr, i);
      }
    }

  } else {

    // Use the work estimate from the last integration cycle.
    result = dataBase.globalWork();

  }

  // We don't allow zeros for work, so if there are any reset them to be the minimum
  // work for a node in the NodeList in question.
  const double globalMax = result.max();
  if (globalMax == 0.0) {
    result = 1.0;
  } else {
    for (auto iNodeList = 0u; iNodeList != result.size(); ++iNodeList) {
      Field<Dimension, Scalar>& worki = *(result[iNodeList]);
      if (worki.max() == 0.0) {
        worki = globalMax;
      } else {
        if (worki.min() == 0.0) {
          double localMin = std::numeric_limits<double>::max();
          for (auto i = 0u; i != worki.numInternalElements(); ++i) {
            if (worki(i) > 0.0) localMin = std::min(localMin, worki(i));
          }
          double globalMin;
          MPI_Allreduce(&localMin, &globalMin, 1, MPI_DOUBLE, MPI_MIN, Communicator::communicator());
          worki.applyMin(globalMin);
        }
      }
    }
  }

  // Output some statistics.
  const Scalar minWeight = result.min();
  const Scalar maxWeight = result.max();
  if (Process::getRank() == 0) cout << "RedistributeNodes::workPerNode: min/max work : "
                                    << minWeight << " "
                                    << maxWeight << endl;

  // Return the result.
  return result;
}

//------------------------------------------------------------------------------
// Pack a vector<DomainNode> into vector<char>.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<char>
RedistributeNodes<Dimension>::
packDomainNodes(const vector<DomainNode<Dimension> >& domainNodes) const {
  vector<char> result;
  const int bufsize = (4*sizeof(int) + Dimension::nDim*sizeof(double))*domainNodes.size();
  result.reserve(bufsize);
  for (typename vector<DomainNode<Dimension> >::const_iterator domainNodeItr = domainNodes.begin();
       domainNodeItr != domainNodes.end();
       ++domainNodeItr) {
    packElement(domainNodeItr->localNodeID, result);
    packElement(domainNodeItr->globalNodeID, result);
    packElement(domainNodeItr->nodeListID, result);
    packElement(domainNodeItr->domainID, result);
    for (int i = 0; i != Dimension::nDim; ++i)
      packElement(domainNodeItr->position(i), result);
  }
  ENSURE((int)result.size() == bufsize);
  return result;
}

//------------------------------------------------------------------------------
// Unpack a vector<char> to a vector<DomainNode>.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<DomainNode<Dimension> >
RedistributeNodes<Dimension>::
unpackDomainNodes(const vector<char>& packedDomainNodes) const {
  const int sizeOfDomainNode = 4*sizeof(int) + Dimension::nDim*sizeof(double);
  const int size = packedDomainNodes.size()/sizeOfDomainNode;
  CHECK(packedDomainNodes.size() - size*sizeOfDomainNode == 0);

  vector<DomainNode<Dimension> > result;
  result.reserve(size);
  vector<char>::const_iterator itr = packedDomainNodes.begin();
  const vector<char>::const_iterator endItr = packedDomainNodes.end();
  while (itr < packedDomainNodes.end()) {
    result.push_back(DomainNode<Dimension>());
    unpackElement(result.back().localNodeID, itr, endItr);
    unpackElement(result.back().globalNodeID, itr, endItr);
    unpackElement(result.back().nodeListID, itr, endItr);
    unpackElement(result.back().domainID, itr, endItr);
    for (int i = 0; i < Dimension::nDim; ++i) 
      unpackElement(result.back().position(i), itr, endItr);
  }
  ENSURE(itr == packedDomainNodes.end());
  return result;
}

//------------------------------------------------------------------------------
// Gather statistics about the work distribution, returning a string with the
// result.
//------------------------------------------------------------------------------
template<typename Dimension>
string
RedistributeNodes<Dimension>::
gatherDomainDistributionStatistics(const FieldList<Dimension, typename Dimension::Scalar>& work) const {

  // Each domain computes it's total work and number of nodes.
  int localNumNodes = 0;
  Scalar localWork = 0.0;
  for (InternalNodeIterator<Dimension> nodeItr = work.internalNodeBegin();
       nodeItr != work.internalNodeEnd();
       ++nodeItr) {
    ++localNumNodes;
    localWork += work(nodeItr);
  }

  // Now gather up some statistics about the distribution.
  const int numProcs = this->numDomains();
  CHECK(numProcs > 0);
  int globalMinNodes, globalMaxNodes, globalAvgNodes;
  MPI_Allreduce(&localNumNodes, &globalMinNodes, 1, MPI_INT, MPI_MIN, Communicator::communicator());
  MPI_Allreduce(&localNumNodes, &globalMaxNodes, 1, MPI_INT, MPI_MAX, Communicator::communicator());
  MPI_Allreduce(&localNumNodes, &globalAvgNodes, 1, MPI_INT, MPI_SUM, Communicator::communicator());
  globalAvgNodes /= numProcs;

  Scalar globalMinWork, globalMaxWork, globalAvgWork;
  MPI_Allreduce(&localWork, &globalMinWork, 1, MPI_DOUBLE, MPI_MIN, Communicator::communicator());
  MPI_Allreduce(&localWork, &globalMaxWork, 1, MPI_DOUBLE, MPI_MAX, Communicator::communicator());
  MPI_Allreduce(&localWork, &globalAvgWork, 1, MPI_DOUBLE, MPI_SUM, Communicator::communicator());
  globalAvgWork /= numProcs;

  // Build a string with the result.
  std::stringstream result;
  result << "    (min, max, avg) nodes per domain: ("
         << globalMinNodes << ", "
         << globalMaxNodes << ", "
         << globalAvgNodes << ")" << endl
         << "    (min, max, avg) work per domain : ("
         << globalMinWork << ", "
         << globalMaxWork << ", "
         << globalAvgWork << ")" << std::endl;
  return result.str();
}

}

