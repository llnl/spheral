//---------------------------------Spheral++----------------------------------//
// globalNodeIDs (contains numGlobalNodes and globalNodeIDs)
//
// Helper methods to assign unique global node IDs to all node in a NodeList,
// serial or parallel.  *Should* compute the same unique ID for each node no
// matter how the parallel domains are carved up.  
//
// Note this is accomplished by sorting based on position, so nodes occupying
// the same position are considered degenerate and may get different IDs for
// different decompositions.
//
// Created by JMO, Fri Oct 29 13:08:33 2004
//----------------------------------------------------------------------------//
#include "NodeList/NodeList.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "Utilities/peanoHilbertOrderIndices.hh"
#include "Utilities/KeyTraits.hh"
#include "Utilities/DBC.hh"
#include "Utilities/range.hh"
#include "Distributed/allReduce.hh"
#include "Distributed/Communicator.hh"

#include <vector>
#include <tuple>

#ifdef SPHERAL_ENABLE_MPI
#include <mpi.h>
#endif

namespace Spheral {

//------------------------------------------------------------------------------
// Return the total global number of nodes in the NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t
numGlobalNodes(const NodeList<Dimension>& nodeList) {
  auto localResult = nodeList.numInternalNodes();
  return allReduce(localResult, SPHERAL_OP_SUM);
}
  
//------------------------------------------------------------------------------
// Return the total global number of nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t
numGlobalNodes(const DataBase<Dimension>& dataBase) {
  return numGlobalNodes<Dimension, typename DataBase<Dimension>::ConstNodeListIterator>(dataBase.nodeListBegin(), dataBase.nodeListEnd());
}

//------------------------------------------------------------------------------
// Return the total global number of nodes in the set of NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension, typename NodeListIterator>
inline
size_t
numGlobalNodes(const NodeListIterator& begin,
               const NodeListIterator& end) {
  size_t result = 0u;
  for (auto nodeListItr = begin; nodeListItr != end; ++nodeListItr) {
    result += numGlobalNodes(**nodeListItr);
  }
  return result;
}
  
//------------------------------------------------------------------------------
// Compute a unique set of global node IDs for the given NodeList, and return
// the set of them on this process.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, size_t>
globalNodeIDs(const NodeList<Dimension>& nodeList) {

  typedef typename KeyTraits::Key Key;

  // Get the local domain ID and number of processors.
  const size_t procID = Process::getRank();
  const size_t numProcs = Process::getTotalNumberOfProcesses();

  // Build keys to sort the nodes by.
  DataBase<Dimension> db;
  db.appendNodeList(const_cast<NodeList<Dimension>&>(nodeList));
  FieldList<Dimension, Key> keys = peanoHilbertOrderIndices(db);

  // Build the local list of node info.
  typedef std::vector<std::tuple<Key, size_t, size_t> > InfoType;
  InfoType nodeInfo;
  auto numLocalNodes = nodeList.numInternalNodes();
  for (auto i = 0u; i < numLocalNodes; ++i) {
    nodeInfo.push_back(std::make_tuple(keys(0, i), i, procID));
  }

  // Reduce the list of node info to processor 0.
  auto nglobal = numLocalNodes;
  CONTRACT_VAR(nglobal);
#ifdef SPHERAL_ENABLE_MPI
  if (procID == 0u) {

    // Process 0 receives and builds the global info.
    for (auto recvDomain = 1u; recvDomain < numProcs; ++recvDomain) {
      MPI_Status status;
      size_t bufsize;
      MPI_Recv(&bufsize, 1, DataTypeTraits<size_t>::MpiDataType(), recvDomain, 10, Communicator::communicator(), &status);
      CHECK(bufsize > 0u);
      std::vector<char> buf(bufsize);
      MPI_Recv(&buf.front(), bufsize, MPI_CHAR, recvDomain, 11, Communicator::communicator(), &status);

      size_t numRecvNodes;
      std::vector<char>::const_iterator bufItr = buf.begin();
      unpackElement(numRecvNodes, bufItr, buf.end());
      Key key;
      size_t id;
      for (auto i = 0u; i < numRecvNodes; ++i) {
        unpackElement(key, bufItr, buf.end());
        unpackElement(id, bufItr, buf.end());
        nodeInfo.emplace_back(key, id, recvDomain);
      }
      nglobal += numRecvNodes;
      CHECK(bufItr == buf.end());
    }
    CHECK(nodeInfo.size() == nglobal);

  } else {
             
    // Send our local info to processor 0.
    std::vector<char> buf;
    packElement(numLocalNodes, buf);
    for (const auto& [key, id, proc]: nodeInfo) {
      packElement(key, buf);
      packElement(id, buf);
    }
    size_t bufsize = buf.size();
    MPI_Send(&bufsize, 1, DataTypeTraits<size_t>::MpiDataType(), 0, 10, Communicator::communicator());
    MPI_Send(&buf.front(), bufsize, MPI_CHAR, 0, 11, Communicator::communicator());
  }

  CHECK((procID == 0u and nodeInfo.size() == nglobal) or nodeInfo.size() == numLocalNodes);
#endif

  // Sort the node info.
  if (nodeInfo.size() > 0) {
    std::sort(nodeInfo.begin(), nodeInfo.end());
    BEGIN_CONTRACT_SCOPE
    for (auto i = 0u; i < nodeInfo.size() - 1u; ++i) {
      CHECK(std::get<0>(nodeInfo[i]) <= std::get<0>(nodeInfo[i + 1]));
    }
    END_CONTRACT_SCOPE
  }

  // Now we can assign consecutive global IDs based on the sorted list.
  std::vector<std::vector<size_t>> globalIDs(numProcs);
  if (procID == 0u) {
    size_t iglobal = 0u;
    for (const auto& [key, localID, recvProc]: nodeInfo) {
      CHECK(recvProc < globalIDs.size());
      if (localID + 1u > globalIDs[recvProc].size()) globalIDs[recvProc].resize(localID + 1u);
      globalIDs[recvProc][localID] = iglobal++;
    }
    CHECK(iglobal == nglobal);
  }
  Barrier(Communicator::communicator());

  // Assign process 0's Ids.
  Field<Dimension, size_t> result("global IDs", nodeList);
  if (procID == 0u) {
    for (auto i = 0u; i < numLocalNodes; ++i) result(i) = globalIDs[0][i];
  }

#ifdef SPHERAL_ENABLE_MPI
  // Farm the ID info back to all processors.
  if (procID == 0u) {

    // Process 0 sends the info.
    for (auto recvProc = 1u; recvProc < numProcs; ++recvProc) {
      size_t numRecvNodes = globalIDs[recvProc].size();
      MPI_Send(&numRecvNodes, 1, DataTypeTraits<size_t>::MpiDataType(), recvProc, 20, Communicator::communicator());
      if (numRecvNodes > 0u) MPI_Send(&(*globalIDs[recvProc].begin()), numRecvNodes, DataTypeTraits<size_t>::MpiDataType(),
                                      recvProc, 21, Communicator::communicator());
    }

  } else {

    // Get our ids from process 0.
    MPI_Status status;
    size_t numRecvNodes;
    MPI_Recv(&numRecvNodes, 1, DataTypeTraits<size_t>::MpiDataType(), 0, 20, Communicator::communicator(), &status);
    CHECK2(numRecvNodes == numLocalNodes, numRecvNodes << " != " << numLocalNodes);
    if (numRecvNodes > 0) MPI_Recv(&(result[0]), numRecvNodes, DataTypeTraits<size_t>::MpiDataType(), 0, 21,
                                   Communicator::communicator(), &status);

  }
#endif

  // We're done.
  return result;
}

//------------------------------------------------------------------------------
// Compute a unique set of global node IDs for all nodes across all NodeLists in
// a DataBase, returning the result as a FieldList<int>.
//------------------------------------------------------------------------------
template<typename Dimension>
// inline
FieldList<Dimension, size_t>
globalNodeIDs(const DataBase<Dimension>& dataBase) {
  return globalNodeIDs<Dimension, typename DataBase<Dimension>::ConstNodeListIterator>(dataBase.nodeListBegin(), dataBase.nodeListEnd());
}

//------------------------------------------------------------------------------
// Compute a unique set of global node IDs for all nodes in the given set of
// NodeLists, returning the result as a FieldList<int>.
//------------------------------------------------------------------------------
template<typename Dimension, typename NodeListIterator>
// inline
FieldList<Dimension, size_t>
globalNodeIDs(const NodeListIterator& begin,
              const NodeListIterator& end) {

  // Prepare the result.
  const size_t numNodeLists = std::distance(begin, end);
  CONTRACT_VAR(numNodeLists);
  FieldList<Dimension, size_t> result(FieldStorageType::CopyFields);
  for (auto* nodeListPtr: range(begin, end)) result.appendNewField("global IDs", *nodeListPtr, size_t(0));
  CHECK(result.numFields() == numNodeLists);

#ifdef SPHERAL_ENABLE_MPI
  // This processors domain id.
  const size_t procID = Process::getRank();
  const size_t numProcs = Process::getTotalNumberOfProcesses();

  // Count up how many nodes are on this domain.
  size_t numDomainNodes = 0u;
  for (auto nodeListItr = begin; nodeListItr != end; ++nodeListItr)
    numDomainNodes += (*nodeListItr)->numInternalNodes();

  // Now loop over each processor, and have each send the cumulative number
  // on to the next.
  size_t beginID = 0;
  for (auto sendProc = 0u; sendProc < numProcs - 1u; ++sendProc) {
    size_t sendProcDomainNodes = numDomainNodes;
    MPI_Bcast(&sendProcDomainNodes, 1, DataTypeTraits<size_t>::MpiDataType(), sendProc, Communicator::communicator());
    if (procID > sendProc) beginID += sendProcDomainNodes;
  }
  const size_t endID = beginID + numDomainNodes;

  // Now assign global IDs to each node on this domain in the DataBase.
  // Loop over each NodeList in the DataBase.
  size_t nodeListID = 0u;
  for (auto nodeListItr = begin; nodeListItr != end; ++nodeListItr, ++nodeListID) {
    const auto& nodeList = **nodeListItr;
    Field<Dimension, size_t>& globalIDs = **result.fieldForNodeList(nodeList);

    // Construct the set of global IDs for this NodeList on this domain.
    CONTRACT_VAR(endID);
    CHECK(endID - beginID >= nodeList.numInternalNodes());
    for (auto i = 0u; i < nodeList.numInternalNodes(); ++i) globalIDs(i) = beginID + i;
    beginID += nodeList.numInternalNodes();
  }
  CHECK(beginID == endID);

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    ENSURE(result.numFields() == numNodeLists);

    // Make sure each global ID is unique.
    const size_t nGlobal = numGlobalNodes<Dimension, NodeListIterator>(begin, end);
    for (auto checkProc = 0u; checkProc < numProcs; ++checkProc) {
      for (auto* fieldPtr: result) {
        size_t n = fieldPtr->nodeListPtr()->numInternalNodes();
        auto fieldBegin = fieldPtr->begin();
        auto fieldEnd = fieldBegin + n;
        MPI_Bcast(&n, 1, DataTypeTraits<size_t>::MpiDataType(), checkProc, Communicator::communicator());
        for (size_t i = 0u; i < n; ++i) {
          size_t id;
          if (procID == checkProc) id = (*fieldPtr)(i);
          MPI_Bcast(&id, 1, DataTypeTraits<size_t>::MpiDataType(), checkProc, Communicator::communicator());
          CONTRACT_VAR(nGlobal);
          ENSURE(i < nGlobal);
          if (procID != checkProc){
            CONTRACT_VAR(fieldEnd);
            ENSURE(find(fieldBegin, fieldEnd, id) == fieldEnd);
          }
        }
      }
    }
  }
  END_CONTRACT_SCOPE

#else

  // Provide a serial version (pretty simple mindedly).
  size_t numCumulativeNodes = 0u;
  for (auto nodeListItr = begin; nodeListItr != end; ++nodeListItr) {
    const auto& nodeList = **nodeListItr;
    Field<Dimension, size_t>& globalIDs = **result.fieldForNodeList(nodeList);
    globalIDs = globalNodeIDs(nodeList);
    for (auto i = 0u; i < globalIDs.numElements(); ++i) globalIDs(i) += numCumulativeNodes;
    numCumulativeNodes += globalIDs.numElements();
  }

#endif

  // That's it.
  return result;
}

}
