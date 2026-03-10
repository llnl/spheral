//---------------------------------Spheral++----------------------------------//
// DistributedBoundary -- Base class for distributed parallel boundary
// conditions, connecting NodeLists across parallel domains.
//
// Created by JMO, Mon Aug 27 16:57:11 PDT 2001
//----------------------------------------------------------------------------//
#include "DistributedBoundary.hh"
#include "Communicator.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "Field/Field.hh"
#include "Utilities/DataTypeTraits.hh"
#include "Utilities/packElement.hh"
#include "Utilities/removeElements.hh"
#include "Utilities/range.hh"
#include "Utilities/DBC.hh"
#include "waitAllWithDeadlockDetection.hh"
#include "Process.hh"

#include <sstream>
#include <list>
#include <algorithm>
using std::vector;
using std::list;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
DistributedBoundary<Dimension>::DistributedBoundary():
  Boundary<Dimension>(),
  mDomainID(-1),
  mNodeListDomainBoundaryNodeMap(),
  mExchangeFields(),
  mMPIFieldTag(0),
  mSendRequests(),
  mRecvRequests(),
#ifdef SPHERAL_ENABLE_MPI_DEADLOCK_DETECTION
  mSendProcIDs(),
  mRecvProcIDs(),
#endif
  mSendBuffers(),
  mRecvBuffers() {
  
  // Get the number of processor and this one's rank.
  mDomainID = Process::getRank();
  CHECK(mDomainID >= 0 && mDomainID < numDomains());

  // Reserve space in the request buffers.
  mSendRequests.reserve(100000);
  mRecvRequests.reserve(100000);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
DistributedBoundary<Dimension>::~DistributedBoundary() {
}

//------------------------------------------------------------------------------
// Test if the given NodeList has any domain boundary information.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
DistributedBoundary<Dimension>::
communicatedNodeList(const NodeList<Dimension>& nodeList) const {
  return mNodeListDomainBoundaryNodeMap.find(const_cast<NodeList<Dimension>*>(&nodeList)) !=
    mNodeListDomainBoundaryNodeMap.end();
}

//------------------------------------------------------------------------------
// Test if the given NodeList is shared with the given domain.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
DistributedBoundary<Dimension>::
nodeListSharedWithDomain(const NodeList<Dimension>& nodeList,
                         int neighborDomainID) const {
  typename NodeListDomainBoundaryNodeMap::const_iterator itr =
    mNodeListDomainBoundaryNodeMap.find(const_cast<NodeList<Dimension>*>(&nodeList));

  if (itr == mNodeListDomainBoundaryNodeMap.end()) return false;
  const DomainBoundaryNodeMap& thpt = itr->second;
  return thpt.find(neighborDomainID) != thpt.end();
}

//------------------------------------------------------------------------------
// Get the domainID <-> DomainBoundaryNodes information for the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
const typename DistributedBoundary<Dimension>::DomainBoundaryNodeMap&
DistributedBoundary<Dimension>::
domainBoundaryNodeMap(const NodeList<Dimension>& nodeList) const {
  typename NodeListDomainBoundaryNodeMap::const_iterator itr =
    mNodeListDomainBoundaryNodeMap.find(const_cast<NodeList<Dimension>*>(&nodeList));
  REQUIRE(itr != mNodeListDomainBoundaryNodeMap.end());
  return itr->second;
}

//------------------------------------------------------------------------------
// Get the DomainBoundaryNodes for the given NodeList and neighbor domain.
//------------------------------------------------------------------------------
template<typename Dimension>
const typename DistributedBoundary<Dimension>::DomainBoundaryNodes&
DistributedBoundary<Dimension>::
domainBoundaryNodes(const NodeList<Dimension>& nodeList,
                    int neighborDomainID) const {
  const DomainBoundaryNodeMap& domainBoundary = domainBoundaryNodeMap(nodeList);
  typename DomainBoundaryNodeMap::const_iterator itr = domainBoundary.find(neighborDomainID);
  CHECK(itr != domainBoundary.end());
  return itr->second;
}

//------------------------------------------------------------------------------
// Extract the set of processes we're sending to and receiving from.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributedBoundary<Dimension>::
communicatedProcs(vector<int>& sendProcs,
                  vector<int>& recvProcs) const {

  // Clear out any preexisting info.
  sendProcs = vector<int>();
  recvProcs = vector<int>();

  // Scan the comm maps and flag any processors we're sending/recving with.
  const size_t numProcs = numDomains();
  vector<int> sendFlags(numProcs, 0);
  vector<int> recvFlags(numProcs, 0);
  for (typename NodeListDomainBoundaryNodeMap::const_iterator itr1 = mNodeListDomainBoundaryNodeMap.begin();
       itr1 != mNodeListDomainBoundaryNodeMap.end();
       ++itr1) {
    for (typename DomainBoundaryNodeMap::const_iterator itr2 = itr1->second.begin();
         itr2 != itr1->second.end();
         ++itr2) {
      const int proc = itr2->first;
      CHECK(proc >= 0 and proc < (int)numProcs);
      const DomainBoundaryNodes& domNodes = itr2->second;
      if (domNodes.sendNodes.size() > 0) sendFlags[proc] = 1;
      if (domNodes.receiveNodes.size() > 0) sendFlags[proc] = 1;
    }
  }

  // Fill in the result.
  for (auto k = 0u; k != numProcs; ++k) {
    if (sendFlags[k] == 1) sendProcs.push_back(k);
    if (recvFlags[k] == 1) recvProcs.push_back(k);
  }
}

// //------------------------------------------------------------------------------
// // Generic method to exchange a Field.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// template<typename DataType>
// void
// DistributedBoundary<Dimension>::
// exchangeField(Field<Dimension, DataType>& field) const {

//   // We use a handy trait class to tell us how many elements there are in the
//   // type we're exchanging.
//   typedef typename DataTypeTraits<DataType>::ElementType ElementType;
//   MPI_Datatype MPI_TYPE = DataTypeTraits<ElementType>::MpiDataType();
//   const int numElementsInType = DataTypeTraits<DataType>::numElements();
//   const NodeList<Dimension>* nodeListPtr = field.nodeListPtr();

//   // Pre-conditions.
//   BEGIN_CONTRACT_SCOPE
//   {
//     const int nProcs = numDomains();
//     const int procID = domainID();
//     for (int checkProc = 0; checkProc != nProcs; ++checkProc) {

//       // First check that everyone agrees about the Field we're working on.
//       string name = field.name();
//       int nameSize = name.size();
//       MPI_Bcast(&nameSize, 1, MPI_INT, checkProc, Communicator::communicator());
//       REQUIRE(nameSize == name.size());
//       MPI_Bcast(&(*name.begin()), nameSize, MPI_CHAR, checkProc, Communicator::communicator());
//       REQUIRE(name == field.name());
//       int elementSize = numElementsInType;
//       MPI_Bcast(&elementSize, 1, MPI_INT, checkProc, Communicator::communicator());
//       REQUIRE(elementSize == numElementsInType);

//       // Determine how many send & receive targets this checkProc thinks it has.
//       int numTargetProcs = 0;
//       const DomainBoundaryNodeMap* domNodeMapPtr = 0;
//       if (communicatedNodeList(*nodeListPtr)) {
//         domNodeMapPtr = &(domainBoundaryNodeMap(*nodeListPtr));
//         for (typename DomainBoundaryNodeMap::const_iterator domainItr = domNodeMapPtr->begin();
//              domainItr != domNodeMapPtr->end();
//              ++domainItr) ++numTargetProcs;
//       }
//       MPI_Bcast(&numTargetProcs, 1, MPI_INT, checkProc, Communicator::communicator());
   
//       // Now check that everyone checkProc is communicating with agrees.
//       if (numTargetProcs > 0) {
//         typename DomainBoundaryNodeMap::const_iterator domainItr;
//         if (procID == checkProc) domainItr = domNodeMapPtr->begin();
//         for (int i = 0; i != numTargetProcs; ++i) {
//           int targetProc = -1;
//           int numSendNodes = 0;
//           int numRecvNodes = 0;
//           if (communicatedNodeList(*nodeListPtr) && procID == checkProc ) {
//             REQUIRE(domainItr != domNodeMapPtr->end());
//             targetProc = domainItr->first;
//             const DomainBoundaryNodes& boundNodes = domainItr->second;
//             numSendNodes = boundNodes.sendNodes.size();
//             numRecvNodes = boundNodes.receiveNodes.size();
//             ++domainItr;
//           }
//           MPI_Bcast(&targetProc, 1, MPI_INT, checkProc, Communicator::communicator());
//           MPI_Bcast(&numSendNodes, 1, MPI_INT, checkProc, Communicator::communicator());
//           MPI_Bcast(&numRecvNodes, 1, MPI_INT, checkProc, Communicator::communicator());
//           REQUIRE(targetProc >= -1 && targetProc < nProcs);
//           REQUIRE(numSendNodes >= 0);
//           REQUIRE(numRecvNodes >= 0);
//           REQUIRE(numSendNodes + numRecvNodes >= 0);

//           // Is everything kosher with the target processor?
//           if (procID == targetProc) {
//             REQUIRE(communicatedNodeList(*nodeListPtr));
//             const DomainBoundaryNodeMap& domBoundNodeMap = domainBoundaryNodeMap(*nodeListPtr);
//             typename DomainBoundaryNodeMap::const_iterator itr = domBoundNodeMap.find(checkProc);
//             REQUIRE(itr != domBoundNodeMap.end());
//             const DomainBoundaryNodes& boundNodes = itr->second;
//             REQUIRE(boundNodes.sendNodes.size() == numRecvNodes);
//             REQUIRE(boundNodes.receiveNodes.size() == numSendNodes);
//           }
//         }
//       }
//     }
//   }
//   END_CONTRACT_SCOPE

//   // Get the map of (domain -> send and receive nodes) for the NodeList
//   // of this Field.
//   list< vector<ElementType> > packedSendValues;
//   list< vector<ElementType> > packedRecvValues;
//   if (communicatedNodeList(*nodeListPtr)) {
//     const DomainBoundaryNodeMap& domBoundNodeMap = domainBoundaryNodeMap(*nodeListPtr);

//     // Post all the receives.
//     vector<MPI_Request> recvRequests;
//     recvRequests.reserve(domBoundNodeMap.size());
//     for (typename DomainBoundaryNodeMap::const_iterator domainItr = domBoundNodeMap.begin();
//          domainItr != domBoundNodeMap.end();
//          ++domainItr) {

//       // Get the set of send/recv nodes for this neighbor domain.
//       const DomainBoundaryNodes& boundNodes = domainItr->second;

//       // If there are any recv nodes for this domain...
//       if (boundNodes.receiveNodes.size() > 0) {

//         // Post a non-blocking receive for this domain.
//         int neighborDomainID = domainItr->first;
//         int numValues = boundNodes.receiveNodes.size()*numElementsInType;
//         packedRecvValues.push_back(vector<ElementType>(numValues));
//         recvRequests.push_back(MPI_Request());
//         vector<ElementType>& recvValues = packedRecvValues.back();
//         MPI_Irecv(&(*recvValues.begin()), numValues, MPI_TYPE,
//                   neighborDomainID, 1, Communicator::communicator(), &(recvRequests.back()));
//       }
//     }
//     CHECK(packedRecvValues.size() == recvRequests.size());

//     // Now send all of our send nodes.
//     vector<MPI_Request> sendRequests;
//     sendRequests.reserve(domBoundNodeMap.size());
//     for (typename DomainBoundaryNodeMap::const_iterator domainItr = domBoundNodeMap.begin();
//          domainItr != domBoundNodeMap.end();
//          ++domainItr) {

//       // Get the set of send/recv nodes for this neighbor domain.
//       const DomainBoundaryNodes& boundNodes = domainItr->second;

//       // If there are any send nodes for this domain...
//       if (boundNodes.sendNodes.size() > 0) {

//         // Do a non-blocking send for this domain.
//         int neighborDomainID = domainItr->first;
//         int numValues = boundNodes.sendNodes.size()*numElementsInType;
//         packedSendValues.push_back(vector<ElementType>(numValues));
//         sendRequests.push_back(MPI_Request());
//      vector<ElementType>& sendValues = packedSendValues.back();
//         packFieldValues(field, boundNodes.sendNodes, sendValues);
//         MPI_Isend(&(*sendValues.begin()), numValues, MPI_TYPE,
//                   neighborDomainID, 1, Communicator::communicator(), &(sendRequests.back()));
//       }
//     }
//     CHECK(packedSendValues.size() == sendRequests.size());

//     // Wait until all of our receives have been satisfied.
//     vector<MPI_Status> recvStatus(recvRequests.size());
//     MPI_Waitall(recvRequests.size(), &(*recvRequests.begin()), &(*recvStatus.begin()));

//     // Unpack encoded Field values into the Field.
//     typename list< vector<ElementType> >::const_iterator recvVectorItr = packedRecvValues.begin();
//     for (typename DomainBoundaryNodeMap::const_iterator domainItr = domBoundNodeMap.begin();
//          domainItr != domBoundNodeMap.end();
//          ++domainItr) {

//       // Get the set of send/recv nodes for this neighbor domain.
//       const DomainBoundaryNodes& boundNodes = domainItr->second;

//       // If there are any recv nodes for this domain...
//       if (boundNodes.receiveNodes.size() > 0) {

//         // The encoded receive values from this domain.
//         CHECK(recvVectorItr != packedRecvValues.end());
//         const vector<ElementType>& recvValues = *recvVectorItr;

//         // Unpack the values from this domain into the Field.
//         unpackFieldValues(field, boundNodes.receiveNodes, recvValues);
//         ++recvVectorItr;
//       }
//     }
//     CHECK(recvVectorItr == packedRecvValues.end());

//     // Wait until all of our sends have been satisfied.
//     vector<MPI_Status> sendStatus(sendRequests.size());
//     MPI_Waitall(sendRequests.size(), &(*sendRequests.begin()), &(*sendStatus.begin()));

//   }

// }

// //------------------------------------------------------------------------------
// // Specialized method to exchange Field<vector<double> >.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void
// DistributedBoundary<Dimension>::
// exchangeField(Field<Dimension, vector<double> >& field) const {

//   // Get the NodeList.
//   const NodeList<Dimension>* nodeListPtr = field.nodeListPtr();

//   // Pre-conditions.
//   BEGIN_CONTRACT_SCOPE
//   {
//     const int nProcs = numDomains();
//     const int procID = domainID();
//     for (int checkProc = 0; checkProc != nProcs; ++checkProc) {

//       // First check that everyone agrees about the Field we're working on.
//       string name = field.name();
//       int nameSize = name.size();
//       MPI_Bcast(&nameSize, 1, MPI_INT, checkProc, Communicator::communicator());
//       REQUIRE(nameSize == name.size());
//       MPI_Bcast(&(*name.begin()), nameSize, MPI_CHAR, checkProc, Communicator::communicator());
//       REQUIRE(name == field.name());

//       // Determine how many send & receive targets this checkProc thinks it has.
//       int numTargetProcs = 0;
//       const DomainBoundaryNodeMap* domNodeMapPtr = 0;
//       if (communicatedNodeList(*nodeListPtr)) {
//         domNodeMapPtr = &(domainBoundaryNodeMap(*nodeListPtr));
//         for (typename DomainBoundaryNodeMap::const_iterator domainItr = domNodeMapPtr->begin();
//              domainItr != domNodeMapPtr->end();
//              ++domainItr) ++numTargetProcs;
//       }
//       MPI_Bcast(&numTargetProcs, 1, MPI_INT, checkProc, Communicator::communicator());
   
//       // Now check that everyone checkProc is communicating with agrees.
//       if (numTargetProcs > 0) {
//         typename DomainBoundaryNodeMap::const_iterator domainItr;
//         if (procID == checkProc) domainItr = domNodeMapPtr->begin();
//         for (int i = 0; i != numTargetProcs; ++i) {
//           int targetProc = -1;
//           int numSendNodes = 0;
//           int numRecvNodes = 0;
//           if (communicatedNodeList(*nodeListPtr) && procID == checkProc ) {
//             REQUIRE(domainItr != domNodeMapPtr->end());
//             targetProc = domainItr->first;
//             const DomainBoundaryNodes& boundNodes = domainItr->second;
//             numSendNodes = boundNodes.sendNodes.size();
//             numRecvNodes = boundNodes.receiveNodes.size();
//             ++domainItr;
//           }
//           MPI_Bcast(&targetProc, 1, MPI_INT, checkProc, Communicator::communicator());
//           MPI_Bcast(&numSendNodes, 1, MPI_INT, checkProc, Communicator::communicator());
//           MPI_Bcast(&numRecvNodes, 1, MPI_INT, checkProc, Communicator::communicator());
//           REQUIRE(targetProc >= -1 && targetProc < nProcs);
//           REQUIRE(numSendNodes >= 0);
//           REQUIRE(numRecvNodes >= 0);
//           REQUIRE(numSendNodes + numRecvNodes >= 0);

//           // Is everything kosher with the target processor?
//           if (procID == targetProc) {
//             REQUIRE(communicatedNodeList(*nodeListPtr));
//             const DomainBoundaryNodeMap& domBoundNodeMap = domainBoundaryNodeMap(*nodeListPtr);
//             typename DomainBoundaryNodeMap::const_iterator itr = domBoundNodeMap.find(checkProc);
//             REQUIRE(itr != domBoundNodeMap.end());
//             const DomainBoundaryNodes& boundNodes = itr->second;
//             REQUIRE(boundNodes.sendNodes.size() == numRecvNodes);
//             REQUIRE(boundNodes.receiveNodes.size() == numSendNodes);
//           }
//         }
//       }
//     }
//   }
//   END_CONTRACT_SCOPE

//   // Get the map of (domain -> send and receive nodes) for the NodeList
//   // of this Field.
//   list< vector<int> > packedNumSendValues;
//   list< vector<int> > packedNumRecvValues;
//   list< vector<double> > packedSendValues;
//   list< vector<double> > packedRecvValues;
//   if (communicatedNodeList(*nodeListPtr)) {
//     const DomainBoundaryNodeMap& domBoundNodeMap = domainBoundaryNodeMap(*nodeListPtr);

//     // Post all of our sends.
//     vector<MPI_Request> sendRequests;
//     sendRequests.reserve(2*domBoundNodeMap.size());
//     for (typename DomainBoundaryNodeMap::const_iterator domainItr = domBoundNodeMap.begin();
//          domainItr != domBoundNodeMap.end();
//          ++domainItr) {

//       // Get the set of send/recv nodes for this neighbor domain.
//       const DomainBoundaryNodes& boundNodes = domainItr->second;

//       // If there are any send nodes for this domain...
//       if (boundNodes.sendNodes.size() > 0) {
//         const int neighborDomainID = domainItr->first;
//         const int numSendNodes = boundNodes.sendNodes.size();

//         // Prepare the send data.
//         packedNumSendValues.push_back(vector<int>(numSendNodes));
//         packedSendValues.push_back(vector<double>());
//         vector<int>& numSendValues = packedNumSendValues.back();
//         vector<double>& sendValues = packedSendValues.back();
//         int totalNumSendValues = 0;
//         for (size_t j = 0; j != numSendNodes; ++j) {
//           CHECK(j < boundNodes.sendNodes.size());
//           CHECK(j < numSendValues.size());
//           const int i = boundNodes.sendNodes[j];
//           numSendValues[j] = field(i).size();
//           copy(field(i).begin(), field(i).end(), back_inserter(sendValues));
//           totalNumSendValues += numSendValues[j];
//         }
//         CHECK(sendValues.size() == totalNumSendValues);

//         // Post a non-blocking send of the number of value per node we're sending.
//         sendRequests.push_back(MPI_Request());
//         MPI_Isend(&(*numSendValues.begin()), numSendNodes, MPI_INT,
//                   neighborDomainID, 1, Communicator::communicator(), &(sendRequests.back()));

//         // Post a non-blocking send of the data.
//         sendRequests.push_back(MPI_Request());
//         MPI_Isend(&(*sendValues.begin()), totalNumSendValues, MPI_DOUBLE, 
//                   neighborDomainID, 2, Communicator::communicator(), &(sendRequests.back()));
//       }
//     }
//     CHECK(packedSendValues.size() == sendRequests.size()/2);
//     CHECK(packedNumSendValues.size() == sendRequests.size()/2);

//     // Post all the receives.
//     vector<MPI_Request> recvRequests;
//     recvRequests.reserve(domBoundNodeMap.size());
//     for (typename DomainBoundaryNodeMap::const_iterator domainItr = domBoundNodeMap.begin();
//          domainItr != domBoundNodeMap.end();
//          ++domainItr) {

//       // Get the set of send/recv nodes for this neighbor domain.
//       const DomainBoundaryNodes& boundNodes = domainItr->second;

//       // If there are any recv nodes for this domain...
//       if (boundNodes.receiveNodes.size() > 0) {
//         const int neighborDomainID = domainItr->first;
//         const int numRecvNodes = boundNodes.receiveNodes.size();

//         // Post a blocking receive for the number of values we're getting
//         // from this domain.
//         MPI_Status recvStatus;
//         packedNumRecvValues.push_back(vector<int>(numRecvNodes));
//         vector<int>& numRecvValues = packedNumRecvValues.back();
//         MPI_Recv(&(*numRecvValues.begin()), numRecvNodes, MPI_INT,
//                   neighborDomainID, 1, Communicator::communicator(), &recvStatus);

//         // Prepare the buffer for the receive values.
//         int totalNumRecvValues = 0;
//         for (typename vector<int>::const_iterator itr = numRecvValues.begin();
//              itr != numRecvValues.end();
//              ++itr) totalNumRecvValues += *itr;
//         packedRecvValues.push_back(vector<double>(totalNumRecvValues));
//         vector<double>& recvValues = packedRecvValues.back();
//         CHECK(recvValues.size() == totalNumRecvValues);

//         // Post a non-blocking receive for the packed values.
//         recvRequests.push_back(MPI_Request());
//         MPI_Irecv(&(*recvValues.begin()), totalNumRecvValues, MPI_DOUBLE,
//                   neighborDomainID, 2, Communicator::communicator(), &(recvRequests.back()));
//       }
//     }
//     CHECK(packedNumRecvValues.size() == packedRecvValues.size());
//     CHECK(packedRecvValues.size() == recvRequests.size());

//     // Wait until all of our receives have been satisfied.
//     vector<MPI_Status> recvStatus(recvRequests.size());
//     MPI_Waitall(recvRequests.size(), &(*recvRequests.begin()), &(*recvStatus.begin()));

//     // Unpack encoded Field values into the Field.
//     typename list< vector<int> >::const_iterator numRecvItr = packedNumRecvValues.begin();
//     typename list< vector<double> >::const_iterator recvVectorItr = packedRecvValues.begin();
//     for (typename DomainBoundaryNodeMap::const_iterator domainItr = domBoundNodeMap.begin();
//          domainItr != domBoundNodeMap.end();
//          ++domainItr) {

//       // Get the set of send/recv nodes for this neighbor domain.
//       const DomainBoundaryNodes& boundNodes = domainItr->second;

//       // If there are any recv nodes for this domain...
//       if (boundNodes.receiveNodes.size() > 0) {
//         const size_t numRecvNodes = boundNodes.receiveNodes.size();

//         // The encoded receive values from this domain.
//         CHECK(numRecvItr != packedNumRecvValues.end());
//         CHECK(recvVectorItr != packedRecvValues.end());
//         const vector<int>& numRecvValues = *numRecvItr;
//         const vector<double>& recvValues = *recvVectorItr;

//         // Unpack the values from this domain into the Field.
//         int offset = 0;
//         for (size_t j = 0; j != numRecvNodes; ++j) {
//           CHECK(j < boundNodes.receiveNodes.size());
//           const int i = boundNodes.receiveNodes[j];
//           CHECK(i >= field.nodeList().firstGhostNode() && i < field.nodeList().numNodes());
//           const int numValues = numRecvValues[j];
//           CHECK(offset + numValues <= recvValues.size());
//           field(i) = vector<double>();
//           copy(recvValues.begin() + offset,
//                recvValues.begin() + offset + numValues,
//                back_inserter(field(i)));
//           offset += numValues;
//         }
//         CHECK(offset == recvValues.size());
//         ++numRecvItr;
//         ++recvVectorItr;
//       }
//     }
//     CHECK(numRecvItr == packedNumRecvValues.end());
//     CHECK(recvVectorItr == packedRecvValues.end());

//     // Wait until all of our sends have been satisfied.
//     vector<MPI_Status> sendStatus(sendRequests.size());
//     MPI_Waitall(sendRequests.size(), &(*sendRequests.begin()), &(*sendStatus.begin()));

//   }

// }

//------------------------------------------------------------------------------
// Begin exchange of a Field with a fixed size DataType.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributedBoundary<Dimension>::
beginExchangeFieldFixedSize(FieldBase<Dimension>& field) const {

  // We use a handy trait class to tell us how many elements there are in the
  // type we're exchanging.
  VERIFY2(field.fixedSizeDataType(), "Assuming we're communicating fixed size types!");
  const int numElementsInType = field.numValsInDataType();
  const NodeList<Dimension>* nodeListPtr = field.nodeListPtr();

  // Grab the number of domains and processor ID.
  const int nProcs = numDomains();
  const int procID = domainID();

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    for (int checkProc = 0; checkProc != nProcs; ++checkProc) {

      // First check that everyone agrees about the Field we're working on.
      string name = field.name();
      int nameSize = name.size();
      MPI_Bcast(&nameSize, 1, MPI_INT, checkProc, Communicator::communicator());
      REQUIRE(nameSize == (int)name.size());
      MPI_Bcast(&(*name.begin()), nameSize, MPI_CHAR, checkProc, Communicator::communicator());
      REQUIRE(name == field.name());
      int elementSize = numElementsInType;
      MPI_Bcast(&elementSize, 1, MPI_INT, checkProc, Communicator::communicator());
      REQUIRE(elementSize == numElementsInType);

      // Determine how many send & receive targets this checkProc thinks it has.
      int numTargetProcs = 0;
      const DomainBoundaryNodeMap* domNodeMapPtr = 0;
      if (communicatedNodeList(*nodeListPtr)) {
        domNodeMapPtr = &(domainBoundaryNodeMap(*nodeListPtr));
        for (typename DomainBoundaryNodeMap::const_iterator domainItr = domNodeMapPtr->begin();
             domainItr != domNodeMapPtr->end();
             ++domainItr) ++numTargetProcs;
      }
      MPI_Bcast(&numTargetProcs, 1, MPI_INT, checkProc, Communicator::communicator());
   
      // Now check that everyone checkProc is communicating with agrees.
      if (numTargetProcs > 0) {
        typename DomainBoundaryNodeMap::const_iterator domainItr;
        if (procID == checkProc) domainItr = domNodeMapPtr->begin();
        for (int i = 0; i != numTargetProcs; ++i) {
          int targetProc = -1;
          int numSendNodes = 0;
          int numRecvNodes = 0;
          if (communicatedNodeList(*nodeListPtr) && procID == checkProc ) {
            REQUIRE(domainItr != domNodeMapPtr->end());
            targetProc = domainItr->first;
            const DomainBoundaryNodes& boundNodes = domainItr->second;
            numSendNodes = boundNodes.sendNodes.size();
            numRecvNodes = boundNodes.receiveNodes.size();
            ++domainItr;
          }
          MPI_Bcast(&targetProc, 1, MPI_INT, checkProc, Communicator::communicator());
          MPI_Bcast(&numSendNodes, 1, MPI_INT, checkProc, Communicator::communicator());
          MPI_Bcast(&numRecvNodes, 1, MPI_INT, checkProc, Communicator::communicator());
          REQUIRE(targetProc >= -1 && targetProc < nProcs);
          REQUIRE(numSendNodes >= 0);
          REQUIRE(numRecvNodes >= 0);
          REQUIRE(numSendNodes + numRecvNodes >= 0);

          // Is everything kosher with the target processor?
          if (procID == targetProc) {
            REQUIRE(communicatedNodeList(*nodeListPtr));
            const DomainBoundaryNodeMap& domBoundNodeMap = domainBoundaryNodeMap(*nodeListPtr);
            typename DomainBoundaryNodeMap::const_iterator itr = domBoundNodeMap.find(checkProc);
            REQUIRE(itr != domBoundNodeMap.end());
            const DomainBoundaryNodes& boundNodes = itr->second;
            CONTRACT_VAR(boundNodes);
            REQUIRE((int)boundNodes.sendNodes.size() == numRecvNodes);
            REQUIRE((int)boundNodes.receiveNodes.size() == numSendNodes);
          }
        }
      }
    }
  }
  END_CONTRACT_SCOPE

  // Establish a unique tag for this field's send/recv ops with neighbor domains.
  ++mMPIFieldTag;

  // We only do work if this is a communicated NodeList on this domain.
  // We also skip if this field is already being exchanged.
  if (communicatedNodeList(*nodeListPtr) and
      mField2SendBuffer.find(&field) == mField2SendBuffer.end()) {

    // Get the map of (domain -> send and receive nodes) for the NodeList
    // of this Field.
    const DomainBoundaryNodeMap& domBoundNodeMap = domainBoundaryNodeMap(*nodeListPtr);

    // Allocate buffers to hold the packed send and receive data.
    mSendBuffers.push_back(list< vector<char> >());
    mRecvBuffers.push_back(list< vector<char> >());
    list< vector<char> >& packedSendValues = mSendBuffers.back();
    list< vector<char> >& packedRecvValues = mRecvBuffers.back();
    mField2SendBuffer[&field] = &packedSendValues;
    mField2RecvBuffer[&field] = &packedRecvValues;

    // Post all the receives.
    for (typename DomainBoundaryNodeMap::const_iterator domainItr = domBoundNodeMap.begin();
         domainItr != domBoundNodeMap.end();
         ++domainItr) {

      // Get the set of send/recv nodes for this neighbor domain.
      const DomainBoundaryNodes& boundNodes = domainItr->second;

      // If there are any recv nodes for this domain...
      if (boundNodes.receiveNodes.size() > 0) {

        // Post a non-blocking receive for this domain.
        const int neighborDomainID = domainItr->first;
        const int bufSize = field.computeCommBufferSize(boundNodes.receiveNodes, neighborDomainID, procID);
        packedRecvValues.push_back(vector<char>(bufSize));
        VERIFY(mRecvRequests.size() < mRecvRequests.capacity() - 1);
        mRecvRequests.push_back(MPI_Request());
        vector<char>& recvValues = packedRecvValues.back();
        // cerr << " --> Recieve buffer for " << field.name() << " from " << neighborDomainID << " : " << mField2RecvBuffer[&field] << endl;
        MPI_Irecv(&(*recvValues.begin()), bufSize, MPI_CHAR,
                  neighborDomainID, mMPIFieldTag, Communicator::communicator(), &(mRecvRequests.back()));

#ifdef SPHERAL_ENABLE_MPI_DEADLOCK_DETECTION
        mRecvProcIDs.push_back(neighborDomainID);
#endif

      }
    }

    // Check the receive sizes.
    BEGIN_CONTRACT_SCOPE
    {
      CHECK2(mRecvBuffers.size() == mField2RecvBuffer.size(), mRecvBuffers.size() << " != " << mField2RecvBuffer.size());
      int totalNumRecvs = 0;
      CONTRACT_VAR(totalNumRecvs);
      for (typename list< list< vector<char> > >::const_iterator itr = mRecvBuffers.begin();
           itr != mRecvBuffers.end();
           ++itr) totalNumRecvs += itr->size();
      CHECK((int)mRecvRequests.size() == totalNumRecvs);
    }
    END_CONTRACT_SCOPE

    // Now send all of our send nodes.
    for (typename DomainBoundaryNodeMap::const_iterator domainItr = domBoundNodeMap.begin();
         domainItr != domBoundNodeMap.end();
         ++domainItr) {

      // Get the set of send/recv nodes for this neighbor domain.
      const DomainBoundaryNodes& boundNodes = domainItr->second;

      // If there are any send nodes for this domain...
      if (boundNodes.sendNodes.size() > 0) {

        // Do a non-blocking send for this domain.
        const int neighborDomainID = domainItr->first;
        VERIFY(mSendRequests.size() < mSendRequests.capacity() - 1);
        mSendRequests.push_back(MPI_Request());
        packedSendValues.push_back(field.packValues(boundNodes.sendNodes));
        vector<char>& sendValues = packedSendValues.back();
        MPI_Isend(&(*sendValues.begin()), sendValues.size(), MPI_CHAR,
                  neighborDomainID, mMPIFieldTag, Communicator::communicator(), &(mSendRequests.back()));

#ifdef SPHERAL_ENABLE_MPI_DEADLOCK_DETECTION
        mSendProcIDs.push_back(neighborDomainID);
#endif

      }
    }

    // Check the send sizes.
    BEGIN_CONTRACT_SCOPE
    {
      CHECK(mSendBuffers.size() == mField2SendBuffer.size());
      int totalNumSends = 0;
      CONTRACT_VAR(totalNumSends);
      for (typename list< list< vector<char> > >::const_iterator itr = mSendBuffers.begin();
           itr != mSendBuffers.end();
           ++itr) totalNumSends += itr->size();
      CHECK((int)mSendRequests.size() == totalNumSends);
    }
    END_CONTRACT_SCOPE
  }
}

//------------------------------------------------------------------------------
// Begin exchange for a Field of variable size DataType.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributedBoundary<Dimension>::
beginExchangeFieldVariableSize(FieldBase<Dimension>& field) const {

  const NodeList<Dimension>* nodeListPtr = field.nodeListPtr();

  // Grab the number of domains and processor ID.
  const int nProcs = numDomains();
  const int procID = domainID();

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    for (int checkProc = 0; checkProc != nProcs; ++checkProc) {

      // First check that everyone agrees about the Field we're working on.
      string name = field.name();
      int nameSize = name.size();
      MPI_Bcast(&nameSize, 1, MPI_INT, checkProc, Communicator::communicator());
      REQUIRE(nameSize == (int)name.size());
      MPI_Bcast(&(*name.begin()), nameSize, MPI_CHAR, checkProc, Communicator::communicator());
      REQUIRE(name == field.name());

      // Determine how many send & receive targets this checkProc thinks it has.
      int numTargetProcs = 0;
      const DomainBoundaryNodeMap* domNodeMapPtr = 0;
      if (communicatedNodeList(*nodeListPtr)) {
        domNodeMapPtr = &(domainBoundaryNodeMap(*nodeListPtr));
        for (typename DomainBoundaryNodeMap::const_iterator domainItr = domNodeMapPtr->begin();
             domainItr != domNodeMapPtr->end();
             ++domainItr) ++numTargetProcs;
      }
      MPI_Bcast(&numTargetProcs, 1, MPI_INT, checkProc, Communicator::communicator());
   
      // Now check that everyone checkProc is communicating with agrees.
      if (numTargetProcs > 0) {
        typename DomainBoundaryNodeMap::const_iterator domainItr;
        if (procID == checkProc) domainItr = domNodeMapPtr->begin();
        for (int i = 0; i != numTargetProcs; ++i) {
          int targetProc = -1;
          int numSendNodes = 0;
          int numRecvNodes = 0;
          if (communicatedNodeList(*nodeListPtr) && procID == checkProc ) {
            REQUIRE(domainItr != domNodeMapPtr->end());
            targetProc = domainItr->first;
            const DomainBoundaryNodes& boundNodes = domainItr->second;
            numSendNodes = boundNodes.sendNodes.size();
            numRecvNodes = boundNodes.receiveNodes.size();
            ++domainItr;
          }
          MPI_Bcast(&targetProc, 1, MPI_INT, checkProc, Communicator::communicator());
          MPI_Bcast(&numSendNodes, 1, MPI_INT, checkProc, Communicator::communicator());
          MPI_Bcast(&numRecvNodes, 1, MPI_INT, checkProc, Communicator::communicator());
          REQUIRE(targetProc >= -1 && targetProc < nProcs);
          REQUIRE(numSendNodes >= 0);
          REQUIRE(numRecvNodes >= 0);
          REQUIRE(numSendNodes + numRecvNodes >= 0);

          // Is everything kosher with the target processor?
          if (procID == targetProc) {
            REQUIRE(communicatedNodeList(*nodeListPtr));
            const DomainBoundaryNodeMap& domBoundNodeMap = domainBoundaryNodeMap(*nodeListPtr);
            typename DomainBoundaryNodeMap::const_iterator itr = domBoundNodeMap.find(checkProc);
            REQUIRE(itr != domBoundNodeMap.end());
            const DomainBoundaryNodes& boundNodes = itr->second;
            CONTRACT_VAR(boundNodes);
            REQUIRE((int)boundNodes.sendNodes.size() == numRecvNodes);
            REQUIRE((int)boundNodes.receiveNodes.size() == numSendNodes);
          }
        }
      }
    }
  }
  END_CONTRACT_SCOPE

  // Establish a unique tag for this field's send/recv ops with neighbor domains.
  ++mMPIFieldTag;

  // We only do work if this is a communicated NodeList on this domain.
  // We also skip if this field is already being exchanged.
  if (communicatedNodeList(*nodeListPtr) and
      mField2SendBuffer.find(&field) == mField2SendBuffer.end()) {

    // Get the map of (domain -> send and receive nodes) for the NodeList
    // of this Field.
    const DomainBoundaryNodeMap& domBoundNodeMap = domainBoundaryNodeMap(*nodeListPtr);

    // Allocate buffers to hold the packed send and receive data.
    mSendBuffers.push_back(list< vector<char> >());
    mRecvBuffers.push_back(list< vector<char> >());
    list< vector<char> >& packedSendValues = mSendBuffers.back();
    list< vector<char> >& packedRecvValues = mRecvBuffers.back();
    mField2SendBuffer[&field] = &packedSendValues;
    mField2RecvBuffer[&field] = &packedRecvValues;

    // We also need to send/recv the buffer sizes since this is not a fixed size DataType.
    const size_t ndoms = domBoundNodeMap.size();
    vector<unsigned> sendBufSizes(ndoms, 0U);
    vector<MPI_Request> sendBufSizeRequests;
    sendBufSizeRequests.reserve(ndoms);

    // Launch all of our sends.
    unsigned idom = 0U;
    for (auto domainItr = domBoundNodeMap.begin();
         domainItr != domBoundNodeMap.end();
         ++domainItr, ++idom) {

      // Get the set of send/recv nodes for this neighbor domain.
      const DomainBoundaryNodes& boundNodes = domainItr->second;

      // If there are any send nodes for this domain...
      if (boundNodes.sendNodes.size() > 0) {

        // Do a non-blocking send for this domain.
        const int neighborDomainID = domainItr->first;
        VERIFY(mSendRequests.size() < mSendRequests.capacity() - 1);
        mSendRequests.push_back(MPI_Request());
        packedSendValues.push_back(field.packValues(boundNodes.sendNodes));
        vector<char>& sendValues = packedSendValues.back();

        // Send the size of the buffer.
        sendBufSizes[idom] = sendValues.size();
        sendBufSizeRequests.push_back(MPI_Request());
        MPI_Isend(&sendBufSizes[idom], 1, MPI_UNSIGNED, neighborDomainID, mMPIFieldTag + 65536, 
                  Communicator::communicator(), &(sendBufSizeRequests.back()));

        // Send the packed field data.
        MPI_Isend(&(*sendValues.begin()), sendValues.size(), MPI_CHAR, neighborDomainID, mMPIFieldTag, 
                  Communicator::communicator(), &(mSendRequests.back()));

#ifdef SPHERAL_ENABLE_MPI_DEADLOCK_DETECTION
        mSendProcIDs.push_back(neighborDomainID);
#endif

      }
    }

    // Check the send sizes.
    BEGIN_CONTRACT_SCOPE;
    {
      CHECK(mSendBuffers.size() == mField2SendBuffer.size());
      int totalNumSends = 0;
      CONTRACT_VAR(totalNumSends);
      for (typename list< list< vector<char> > >::const_iterator itr = mSendBuffers.begin();
           itr != mSendBuffers.end();
           ++itr) totalNumSends += itr->size();
      CHECK((int)mSendRequests.size() == totalNumSends);
      CHECK(sendBufSizeRequests.size() <= ndoms);
    }
    END_CONTRACT_SCOPE

    // Post all the receives.
    for (auto domainItr = domBoundNodeMap.begin();
         domainItr != domBoundNodeMap.end();
         ++domainItr) {

      // Get the set of send/recv nodes for this neighbor domain.
      const int neighborDomainID = domainItr->first;
      const DomainBoundaryNodes& boundNodes = domainItr->second;

      // If there are any recv nodes for this domain...
      if (boundNodes.receiveNodes.size() > 0) {

        // Get the incoming buffer size.
        unsigned bufSize;
        MPI_Status status;
        MPI_Recv(&bufSize, 1, MPI_UNSIGNED, neighborDomainID, mMPIFieldTag + 65536,
                 Communicator::communicator(), &status);

        // Post a non-blocking receive for this domain.
        packedRecvValues.push_back(vector<char>(bufSize));
        VERIFY(mRecvRequests.size() < mRecvRequests.capacity() - 1);
        mRecvRequests.push_back(MPI_Request());
        vector<char>& recvValues = packedRecvValues.back();
        MPI_Irecv(&(*recvValues.begin()), bufSize, MPI_CHAR,
                  neighborDomainID, mMPIFieldTag, Communicator::communicator(), &(mRecvRequests.back()));

#ifdef SPHERAL_ENABLE_MPI_DEADLOCK_DETECTION
        mRecvProcIDs.push_back(neighborDomainID);
#endif

      }
    }

    // Wait until all our send buffer sizes have been received.
    if (not sendBufSizeRequests.empty()) {
      vector<MPI_Status> status(sendBufSizeRequests.size());
      MPI_Waitall(sendBufSizeRequests.size(), &sendBufSizeRequests[0], &status[0]);
    }

    // Check the receive sizes.
    BEGIN_CONTRACT_SCOPE
    {
      CHECK2(mRecvBuffers.size() == mField2RecvBuffer.size(), mRecvBuffers.size() << " != " << mField2RecvBuffer.size());
      int totalNumRecvs = 0;
      CONTRACT_VAR(totalNumRecvs);
      for (typename list< list< vector<char> > >::const_iterator itr = mRecvBuffers.begin();
           itr != mRecvBuffers.end();
           ++itr) totalNumRecvs += itr->size();
      CHECK((int)mRecvRequests.size() == totalNumRecvs);
    }
    END_CONTRACT_SCOPE
  }
}

//------------------------------------------------------------------------------
// Override the required setGhostNodes for NodeLists as a no-op.  We prefer
// the DataBase setGhostNodes for DistributedBoundaries, to facilitate more
// efficient parallel exchanges.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributedBoundary<Dimension>::
setGhostNodes(NodeList<Dimension>&) {
  VERIFY(0);
}

//------------------------------------------------------------------------------
// Set the ghost nodes positions and H tensors, but don't recompute the set
// of ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributedBoundary<Dimension>::
updateGhostNodes(NodeList<Dimension>& nodeList) {

  // Exchange the positions and H fields for this NodeList.
  Field<Dimension, Vector>& positions = nodeList.positions();
  Field<Dimension, SymTensor>& Hfield = nodeList.Hfield();
  applyGhostBoundary(positions);
  applyGhostBoundary(Hfield);
}

//------------------------------------------------------------------------------
// The required method to select nodes in violation of this boundary.  A no-op
// for Distributed boundaries.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributedBoundary<Dimension>::
setViolationNodes(NodeList<Dimension>&) {
}

template<typename Dimension>
void
DistributedBoundary<Dimension>::
updateViolationNodes(NodeList<Dimension>&) {
}

//------------------------------------------------------------------------------
// Provide the required ghost node Boundary conditions for Fields.  Just call
// the templated exchangeField method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributedBoundary<Dimension>::
applyGhostBoundary(FieldBase<Dimension>& field) const {
  if (field.fixedSizeDataType()) {
    // cerr << " -->    FIXED SIZE: " << field.name() << endl;
    beginExchangeFieldFixedSize(field);
  } else {
    // cerr << " --> VARIABLE SIZE: " << field.name() << endl;
    beginExchangeFieldVariableSize(field);
  }
  mExchangeFields.push_back(&field);
}

//------------------------------------------------------------------------------
// Override the base method to cull out inactive ghost nodes based on a
// FieldList of flags.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributedBoundary<Dimension>::
cullGhostNodes(const FieldList<Dimension, size_t>& flagSet,
               FieldList<Dimension, size_t>& old2newIndexMap,
               vector<size_t>& numNodesRemoved) {

  typedef typename Boundary<Dimension>::BoundaryNodes BoundaryNodes;
  const int procID = domainID();

  // typedef NodeListRegistrar<Dimension> Registrar;
  // Registrar& registrar = Registrar::instance();
  // CONTRACT_VAR(registrar);
  // REQUIRE((int)numNodesRemoved.size() == registrar.numNodeLists());

  const vector<size_t> numNodesRemovedPreviously(numNodesRemoved);

  // Count how many messages we're going to be communicating.
  // Note we are reversing the normal communication directions here, since the
  // receiving domain is the one that knows which ghost nodes are actually 
  // needed.  This leads to some confusing uses of "send" and "recv" below.  The
  // rule is the data structures (like mNodeListDomainBoundaryNodeMap) are 
  // backwards, while local variables will be declared as send or recv properly
  // for use in this routine.
  size_t numSendMsgs = 0;
  size_t numRecvMsgs = 0;
  for (auto* fPtr: flagSet) {
    const auto& nodeList = fPtr->nodeList();
    if (communicatedNodeList(nodeList)) {
      const auto& domainBoundaryNodes = this->accessDomainBoundaryNodeMap(nodeList);
      for (auto& [proc, domainNodes]: domainBoundaryNodes) {
        if (domainNodes.receiveNodes.size() > 0) ++numSendMsgs;
        if (domainNodes.sendNodes.size() > 0) ++numRecvMsgs;
      }
    }
  }

  // Post receives for flags from the domains that will be sending to us.
  list<vector<int> > recvBuffers;
  vector<MPI_Request> recvRequests;
  recvRequests.reserve(numRecvMsgs);
  int nodeListOff = 0;
  for (auto flagSetItr = flagSet.begin();
       flagSetItr != flagSet.end();
       ++flagSetItr, ++nodeListOff) {
    const NodeList<Dimension>& nodeList = (**flagSetItr).nodeList();
    if (communicatedNodeList(nodeList)) {
      const DomainBoundaryNodeMap& domainBoundaryNodes = this->accessDomainBoundaryNodeMap(nodeList);
      for (auto& [neighborDomain, domainNodes]: domainBoundaryNodes) {
        if (domainNodes.sendNodes.size() > 0) {
          recvBuffers.push_back(vector<int>(domainNodes.sendNodes.size(), -1));
          recvRequests.push_back(MPI_Request());
          vector<int>& buf = recvBuffers.back();
          const int tag = 1000*neighborDomain + nodeListOff;  // This should be safe so long as we don't have more than 1000 NodeLists.
          MPI_Irecv(&(*buf.begin()), domainNodes.sendNodes.size(), MPI_INT, neighborDomain,
                    tag, Communicator::communicator(), &(recvRequests.back()));
        }
      }
    }
  }
  CHECK(recvBuffers.size() == numRecvMsgs);
  CHECK(recvRequests.size() == numRecvMsgs);

  // Iterate over the domain node sets, remove receive nodes, and send our flags to the neighbors.
  list<vector<int> > sendBuffers;
  vector<MPI_Request> sendRequests;
  sendRequests.reserve(numSendMsgs);
  nodeListOff = 0;
  for (auto flagSetItr = flagSet.begin();
       flagSetItr != flagSet.end();
       ++flagSetItr, ++nodeListOff) {
    const NodeList<Dimension>& nodeList = (**flagSetItr).nodeList();
    if (communicatedNodeList(nodeList)) {
      DomainBoundaryNodeMap& domainBoundaryNodes = this->accessDomainBoundaryNodeMap(nodeList);

      // Grab the flags for this NodeList.
      CHECK(flagSet.haveNodeList(nodeList));
      const auto& flags = **flagSet.fieldForNodeList(nodeList);

      // Go over each domain this NodeList is communicating with.
      const BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(const_cast<NodeList<Dimension>&>(nodeList));
      const size_t myFirstGhostNode = (boundaryNodes.ghostNodes.size() > 0 ? boundaryNodes.ghostNodes[0] : nodeList.numNodes()) - numNodesRemoved[nodeListOff];
      size_t newGhostIndex = myFirstGhostNode;
      for (auto& [neighborDomain, domainNodes]: domainBoundaryNodes) {

        // We only care if we're receiving from this domain.
        if (domainNodes.receiveNodes.size() > 0) {

          // Prepare the send buffer for this domain/NodeList.
          sendBuffers.push_back(vector<int>());
          vector<int>& buf = sendBuffers.back();
          buf.reserve(domainNodes.receiveNodes.size());
          for (auto j: domainNodes.receiveNodes) buf.push_back(flags(j));

          // Send the flags for the nodes we're receiving from this domain.
          CHECK(buf.size() == domainNodes.receiveNodes.size());
          const int tag = 1000*procID + nodeListOff;  // This should be safe so long as we don't have more than 1000 NodeLists.
          sendRequests.push_back(MPI_Request());
          MPI_Isend(&(*buf.begin()), buf.size(), MPI_INT, neighborDomain, tag, Communicator::communicator(), &(sendRequests.back()));

          // Cull the receive nodes.
          {
            vector<size_t> newReceiveNodes;
            for (size_t k = 0; k != domainNodes.receiveNodes.size(); ++k) {
              if (flags(domainNodes.receiveNodes[k]) == 1) {
                newReceiveNodes.push_back(newGhostIndex);
                old2newIndexMap(nodeListOff, domainNodes.receiveNodes[k]) = newGhostIndex;
                ++newGhostIndex;
              } else {
                ++numNodesRemoved[nodeListOff];
              }
            }
            domainNodes.receiveNodes = newReceiveNodes;
          }
        }
      }
    }
  }
  CHECK(sendBuffers.size() == numSendMsgs);
  CHECK(sendRequests.size() == numSendMsgs);

  // Wait until we have all the info from our neighbors.
  if (numRecvMsgs > 0) {
    vector<MPI_Status> recvStatus(numRecvMsgs);
    MPI_Waitall(numRecvMsgs, &(*recvRequests.begin()), &(*recvStatus.begin()));
  }

  // Now go through and remove any send nodes that are no longer needed.
  list<vector<int> >::const_iterator bufItr = recvBuffers.begin();
  nodeListOff = 0;
  for (auto flagSetItr = flagSet.begin();
       flagSetItr != flagSet.end();
       ++flagSetItr, ++nodeListOff) {
    const NodeList<Dimension>& nodeList = (**flagSetItr).nodeList();
    if (communicatedNodeList(nodeList)) {
      const BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(const_cast<NodeList<Dimension>&>(nodeList));
      const size_t myFirstGhostNode = (boundaryNodes.ghostNodes.size() > 0 ? boundaryNodes.ghostNodes[0] : nodeList.numNodes()) - numNodesRemovedPreviously[nodeListOff];
      DomainBoundaryNodeMap& domainBoundaryNodes = this->accessDomainBoundaryNodeMap(nodeList);
      for (auto& [neighborDomain, domainNodes]: domainBoundaryNodes) {

        // If there are nodes we're sending to this domain, then we should have some flags
        // waiting for us.
        if (domainNodes.sendNodes.size() > 0) {
          CHECK(bufItr != recvBuffers.end());
          const vector<int>& flags = *bufItr;
          ++bufItr;
        
          // Remove the send nodes we don't need.
          {
            CHECK(flags.size() == domainNodes.sendNodes.size());
            vector<size_t> newSendNodes;
            for (size_t k = 0; k != domainNodes.sendNodes.size(); ++k) {
              if (flags[k] == 1) {
                CONTRACT_VAR(myFirstGhostNode);
                CHECK(old2newIndexMap(nodeListOff, domainNodes.sendNodes[k]) < myFirstGhostNode);
                newSendNodes.push_back(old2newIndexMap(nodeListOff, domainNodes.sendNodes[k]));
              }
            }
            domainNodes.sendNodes = newSendNodes;
          }
        }
      }
    }
  }

  // Update the control and ghost nodes.
  this->setControlAndGhostNodes();

  // Wait for all sends to complete.
  if (numSendMsgs > 0) {
    vector<MPI_Status> sendStatus(numSendMsgs);
    MPI_Waitall(numSendMsgs, &(*sendRequests.begin()), &(*sendStatus.begin()));
  }
}

//------------------------------------------------------------------------------
// Finalize the ghost boundary condition.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributedBoundary<Dimension>::
finalizeGhostBoundary() const {
  const_cast<DistributedBoundary<Dimension>*>(this)->finalizeExchanges();
}

//------------------------------------------------------------------------------
// Unpack the encoded field buffer to the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributedBoundary<Dimension>::
unpackField(FieldBase<Dimension>& field,
            const list< vector<char> >& packedValues) const {

  // Get the domain mappings for the NodeList of this Field.
  const NodeList<Dimension>& nodeList = field.nodeList();
  const DomainBoundaryNodeMap& domBoundNodeMap = domainBoundaryNodeMap(nodeList);
  CHECK(packedValues.size() <= domBoundNodeMap.size());

  // Loop over receive domains.
  typename list< vector<char> >::const_iterator bufferItr = packedValues.begin();
  for (typename DomainBoundaryNodeMap::const_iterator domainItr = domBoundNodeMap.begin();
       domainItr != domBoundNodeMap.end();
       ++domainItr) {

    // Get the set of send/recv nodes for this neighbor domain.
    const DomainBoundaryNodes& boundNodes = domainItr->second;

    // If there are any receive nodes for this domain...
    if (boundNodes.receiveNodes.size() > 0) {
      CHECK(bufferItr != packedValues.end());
      const vector<char>& recvValues = *bufferItr;
      // cerr << " --> UNPACKING " << field.name() << " from " << domainItr->first << " in buffer " << &recvValues << endl;

      // Unpack this domains values...
      field.unpackValues(boundNodes.receiveNodes, recvValues);
      ++bufferItr;

    }
  }
  ENSURE(bufferItr == packedValues.end());
}

//------------------------------------------------------------------------------
// Finalize the exchanges, which means actually perform them sequentially
// and clear out the pending exchange field buffers.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributedBoundary<Dimension>::finalizeExchanges() {

  BEGIN_CONTRACT_SCOPE
  {
    // Make sure everyone has the same number of exchange fields.
    int nFields = mExchangeFields.size();
    MPI_Bcast(&nFields, 1, MPI_INT, 0, Communicator::communicator());
    REQUIRE(nFields == (int)mExchangeFields.size())
    REQUIRE((int)mSendBuffers.size() <= nFields);
    REQUIRE((int)mRecvBuffers.size() <= nFields);
    REQUIRE((int)mField2SendBuffer.size() <= nFields);
    REQUIRE((int)mField2RecvBuffer.size() <= nFields);

    // Count the numbers of send and receive buffers and requests.
    int numSendBuffers = 0;
    int numRecvBuffers = 0;
    CONTRACT_VAR(numSendBuffers);
    CONTRACT_VAR(numRecvBuffers);
    for (typename list< list< vector<char> > >::const_iterator itr = mSendBuffers.begin();
         itr != mSendBuffers.end();
         ++itr) numSendBuffers += itr->size();
    for (typename list< list< vector<char> > >::const_iterator itr = mRecvBuffers.begin();
         itr != mRecvBuffers.end();
         ++itr) numRecvBuffers += itr->size();
    REQUIRE((int)mSendRequests.size() == numSendBuffers);
    REQUIRE((int)mRecvRequests.size() == numRecvBuffers);
  }
  END_CONTRACT_SCOPE

#ifdef SPHERAL_ENABLE_MPI_DEADLOCK_DETECTION
  vector<int> dummyInts;
  vector<MPI_Request> dummyRequests;
  vector<MPI_Status> dummyStatus;
#endif

  // Do we have any data we're waiting to receive?
  if (mRecvRequests.size() > 0) {

    // Wait until all of our receives have been satisfied.
    vector<MPI_Status> recvStatus(mRecvRequests.size());
#ifdef SPHERAL_ENABLE_MPI_DEADLOCK_DETECTION
    waitallWithDeadlockDetection("DistributedBoundary::finalizeExchanges -- RECEIVE waitall",
                                 dummyInts, mRecvProcIDs, dummyRequests, mRecvRequests, dummyStatus, recvStatus, Communicator::communicator());
#else
    MPI_Waitall(mRecvRequests.size(), &(*mRecvRequests.begin()), &(*recvStatus.begin()));
#endif

    // Unpack all field values.
    for (typename vector<FieldBase<Dimension>*>::const_iterator fieldItr = mExchangeFields.begin();
         fieldItr != mExchangeFields.end();
         ++fieldItr) {
      if (mField2RecvBuffer.find(&(**fieldItr)) != mField2RecvBuffer.end()) unpackField(**fieldItr, *mField2RecvBuffer[&(**fieldItr)]);
    }
  }

  // Do we have any data we're waiting to send?
  if (mSendRequests.size() > 0) {

    // Wait until all of our sends have been satisfied.
    vector<MPI_Status> sendStatus(mSendRequests.size());
#ifdef SPHERAL_ENABLE_MPI_DEADLOCK_DETECTION
    waitallWithDeadlockDetection("DistributedBoundary::finalizeExchanges -- RECEIVE waitall",
                                 mSendProcIDs, dummyInts, mSendRequests, dummyRequests, sendStatus, dummyStatus, Communicator::communicator());
#else
    MPI_Waitall(mSendRequests.size(), &(*mSendRequests.begin()), &(*sendStatus.begin()));
#endif

  }

  // Clear out the pending exchange fields.
  mExchangeFields.clear();

  // Reset the buffers.
  mMPIFieldTag = 0;
  mSendRequests = vector<MPI_Request>();
  mRecvRequests = vector<MPI_Request>();
#ifdef SPHERAL_ENABLE_MPI_DEADLOCK_DETECTION
  mSendProcIDs = vector<int>();
  mRecvProcIDs = vector<int>();
#endif
  mSendRequests.reserve(100000);
  mRecvRequests.reserve(100000);
  mSendBuffers = CommBufferSet();
  mRecvBuffers = CommBufferSet();
  mField2SendBuffer = Field2BufferType();
  mField2RecvBuffer = Field2BufferType();

  // Post-conditions.
  ENSURE(mExchangeFields.size() == 0);
  ENSURE(mMPIFieldTag == 0);
  ENSURE(mSendRequests.size() == 0);
  ENSURE(mRecvRequests.size() == 0);
  ENSURE(mSendRequests.capacity() >= 100000);
  ENSURE(mRecvRequests.capacity() >= 100000);
  ENSURE(mSendBuffers.size() == 0);
  ENSURE(mRecvBuffers.size() == 0);
  ENSURE(mField2SendBuffer.size() == 0);
  ENSURE(mField2RecvBuffer.size() == 0);
}

//------------------------------------------------------------------------------
// Update the control and ghost node lists of the base class.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributedBoundary<Dimension>::
setControlAndGhostNodes() {

  const auto& registrar = NodeListRegistrar<Dimension>::instance();
  for (auto* nodeListPtr: range(registrar.begin(), registrar.end())) {

    // Add an entry for this NodeList.
    this->addNodeList(*nodeListPtr);
    auto& boundNodes = this->accessBoundaryNodes(*nodeListPtr);
    auto& controlNodes = boundNodes.controlNodes;
    auto& ghostNodes = boundNodes.ghostNodes;
    controlNodes.clear();
    ghostNodes.clear();

    // Loop over all the domains this NodeList communicates with.
    if (communicatedNodeList(*nodeListPtr)) {
      const auto& domainBoundaryNodes = domainBoundaryNodeMap(*nodeListPtr);
      for (auto& [neighborDomain, domainNodes]: domainBoundaryNodes) {

        // Add the send nodes for this domain to the control node list.
        const vector<size_t>& sendNodes = domainNodes.sendNodes;
        CHECK(controlNodes.size() + sendNodes.size() < 10000000);
        controlNodes.reserve(controlNodes.size() + sendNodes.size());
        for (auto j: sendNodes) controlNodes.push_back(j);

        // Add the receive nodes to the ghost node list.
        const vector<size_t>& recvNodes = domainNodes.receiveNodes;
        CHECK(ghostNodes.size() + sendNodes.size() < 10000000);
        ghostNodes.reserve(ghostNodes.size() + sendNodes.size());
        for (auto j: recvNodes) ghostNodes.push_back(j);
      }      
    }
  }
}

//------------------------------------------------------------------------------
// Read/Write access to the domainID <-> DomainBoundaryNodes information for
// the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
typename DistributedBoundary<Dimension>::DomainBoundaryNodeMap&
DistributedBoundary<Dimension>::
accessDomainBoundaryNodeMap(const NodeList<Dimension>& nodeList) {
  return const_cast<DomainBoundaryNodeMap&>(domainBoundaryNodeMap(nodeList));
}

//------------------------------------------------------------------------------
// Read/Write access to the DomainBoundaryNodes for the given NodeList and
// neighbor domain.
//------------------------------------------------------------------------------
template<typename Dimension>
typename DistributedBoundary<Dimension>::DomainBoundaryNodes&
DistributedBoundary<Dimension>::
accessDomainBoundaryNodes(const NodeList<Dimension>& nodeList,
                          int neighborDomainID) {
  return const_cast<DomainBoundaryNodes&>(domainBoundaryNodes(nodeList,
                                                              neighborDomainID));
}

//------------------------------------------------------------------------------
// Convenience method to return an iterator to the DomainBoundaryNodes for the
// given NodeList and domain pair.  If there isn't an entry for this pair already,
// this method creates one.
//------------------------------------------------------------------------------
template<typename Dimension>
typename DistributedBoundary<Dimension>::DomainBoundaryNodes&
DistributedBoundary<Dimension>::
openDomainBoundaryNodes(NodeList<Dimension>* nodeListPtr,
                        const int domainID) {

  // Is there an entry for this NodeList yet?
  typename NodeListDomainBoundaryNodeMap::iterator nodeListItr = mNodeListDomainBoundaryNodeMap.find(nodeListPtr);

  // If not, create one.
  if (nodeListItr == mNodeListDomainBoundaryNodeMap.end()) {
    mNodeListDomainBoundaryNodeMap[nodeListPtr] = DomainBoundaryNodeMap();
    nodeListItr = mNodeListDomainBoundaryNodeMap.find(nodeListPtr);
    CHECK(nodeListItr != mNodeListDomainBoundaryNodeMap.end());
  }

  // Now is there an entry for the given domain?
  typename DomainBoundaryNodeMap::iterator domainItr = nodeListItr->second.find(domainID);

  // If not, create one.
  if (domainItr == nodeListItr->second.end()) {
    nodeListItr->second[domainID] = DomainBoundaryNodes();
    domainItr = nodeListItr->second.find(domainID);
    CHECK(domainItr != nodeListItr->second.end());
  }

  return domainItr->second;
}

//------------------------------------------------------------------------------
// Inverse of above -- remove the DomainBoundaryNodes for a NodeList/procID pair.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributedBoundary<Dimension>::
removeDomainBoundaryNodes(NodeList<Dimension>* nodeListPtr,
                          const int domainID) {

  typename NodeListDomainBoundaryNodeMap::iterator nodeListItr = mNodeListDomainBoundaryNodeMap.find(nodeListPtr);

  // Is there an entry for this NodeList yet?
  if (nodeListItr != mNodeListDomainBoundaryNodeMap.end()) {

    // Is this domain present?
    typename DomainBoundaryNodeMap::iterator domainItr = nodeListItr->second.find(domainID);
    if (domainItr != nodeListItr->second.end()) nodeListItr->second.erase(domainItr);

    // Are there any domains left for this NodeList?
    if (nodeListItr->second.size() == 0) mNodeListDomainBoundaryNodeMap.erase(nodeListItr);
  }
}

//------------------------------------------------------------------------------
// Clear out any NodeList information that is currently present.
// (Override of the Boundary method.)
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributedBoundary<Dimension>::reset(const DataBase<Dimension>& dataBase) {

  // Call the ancestor method.
  Boundary<Dimension>::reset(dataBase);

  // Clear our own internal data.
  for (typename DataBase<Dimension>::ConstNodeListIterator iter = 
    dataBase.nodeListBegin(); iter != dataBase.nodeListEnd(); ++iter) {
    mNodeListDomainBoundaryNodeMap.erase(*iter);
  } // end for
}

//------------------------------------------------------------------------------
// Distribute the send nodes, filling in every processors receive node info.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributedBoundary<Dimension>::
buildReceiveAndGhostNodes(const DataBase<Dimension>& dataBase) {

  // This processor's ID.
  int procID = this->domainID();
  int numProcs = this->numDomains();

  const int numNodeLists = dataBase.numNodeLists();
  vector<int> firstNewGhostNode(numNodeLists);

  // This if prevents a dereference of a zero length vector &recvStatus.front()
  if (numProcs > 1) {
    CHECK(procID < numProcs);

    // Reserve space for the number of nodes we'll be getting from each domain.
    vector<vector<int> > numSendNodes(numProcs), numRecvNodes(numProcs);
    for (int neighborProc = 0; neighborProc != numProcs; ++neighborProc) {
      numSendNodes[neighborProc].resize(size_t(numNodeLists), 0);
      numRecvNodes[neighborProc].resize(size_t(numNodeLists), 0);
    }

    // Post receives for how many nodes we'll be getting from each domain per
    // NodeList.
    vector<MPI_Request> recvRequests;
    recvRequests.reserve(numProcs - 1);
    for (int neighborProc = 0; neighborProc != numProcs; ++neighborProc) {
      if (neighborProc != procID) {
        recvRequests.push_back(MPI_Request());
        MPI_Irecv(&numRecvNodes[neighborProc].front(), numNodeLists, MPI_INT, neighborProc, 1001, Communicator::communicator(), &(recvRequests.back()));
      }
    }
    CHECK((int)recvRequests.size() == numProcs - 1);

    // Determine how many nodes per NodeList we're sending to each domain.
    for (int neighborProc = 0; neighborProc != numProcs; ++neighborProc) {
      if (neighborProc != procID) {
        int nodeListID = 0;
        for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = dataBase.nodeListBegin();
             nodeListItr != dataBase.nodeListEnd();
             ++nodeListItr, ++nodeListID) {
          if (this->nodeListSharedWithDomain(**nodeListItr, neighborProc)) {
            const vector<size_t>& sendNodes = this->accessDomainBoundaryNodes(**nodeListItr, neighborProc).sendNodes;
            CHECK(sendNodes.size() > 0);
            numSendNodes[neighborProc][nodeListID] = sendNodes.size();
          }
        }
      }
    }

    // Send everyone the sizes.
    for (int neighborProc = 0; neighborProc != numProcs; ++neighborProc) {
      if (neighborProc != procID) {
        MPI_Send(&numSendNodes[neighborProc].front(), numNodeLists, MPI_INT, neighborProc, 1001, Communicator::communicator());
      }
    }

    // Wait until our receives are satisfied.
    // This if prevents a dereference of a zero length vector &recvStatus.front()
    if (numProcs > 1) {
      const int numRecv = numProcs - 1;
      vector<MPI_Status> recvStatus(numRecv);
      MPI_Waitall(numRecv, &recvRequests.front(), &recvStatus.front());
    }

    // Count up the total number of new nodes we'll require for each NodeList.
    vector<int> numNewGhostNodes(size_t(numNodeLists), 0);
    for (int neighborProc = 0; neighborProc != numProcs; ++neighborProc) {
      for (int nodeListID = 0; nodeListID != numNodeLists; ++nodeListID) {
        numNewGhostNodes[nodeListID] += numRecvNodes[neighborProc][nodeListID];
      }
    }

    // Allocate the new ghost nodes on each NodeList in one shot.
    {
      int nodeListID = 0;
      for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = dataBase.nodeListBegin();
           nodeListItr != dataBase.nodeListEnd();
           ++nodeListItr, ++nodeListID) {
        firstNewGhostNode[nodeListID] = (**nodeListItr).numNodes();
        if (numNewGhostNodes[nodeListID] > 0) {
          const int currentNumGhostNodes = (**nodeListItr).numGhostNodes();
          (**nodeListItr).numGhostNodes(currentNumGhostNodes + numNewGhostNodes[nodeListID]);
        }
      }
    }

    // Associate a slice of the new ghost nodes with each communicating domain, filling in the 
    // the receive indicies appropriately.
    for (int neighborProc = 0; neighborProc != numProcs; ++neighborProc) {
      if (neighborProc != procID) {
        int nodeListID = 0;
        for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = dataBase.nodeListBegin();
             nodeListItr != dataBase.nodeListEnd();
             ++nodeListItr, ++nodeListID) {
          if (numRecvNodes[neighborProc][nodeListID] > 0) {
            DomainBoundaryNodes& domainNodes = openDomainBoundaryNodes(&(**nodeListItr), neighborProc);
            vector<size_t>& recvNodes = domainNodes.receiveNodes;
            recvNodes.reserve(numRecvNodes[neighborProc][nodeListID]);
            for (int i = firstNewGhostNode[nodeListID]; 
                 i != firstNewGhostNode[nodeListID] + numRecvNodes[neighborProc][nodeListID];
                 ++i) recvNodes.push_back(i);
            CHECK((int)recvNodes.size() == numRecvNodes[neighborProc][nodeListID]);
            firstNewGhostNode[nodeListID] += numRecvNodes[neighborProc][nodeListID];
          }
        }
      }
    }
  }

  // Fill out the Boundary list of control and ghost nodes with the send and receive
  // nodes on this domain as a courtesy.
  // At the moment this information won't be used for this boundary condition.
  this->setControlAndGhostNodes();

  BEGIN_CONTRACT_SCOPE
  {
    // Ensure we wound up with the correct slices of ghost node indicies.
    if (numProcs >1) {
      int nodeListID = 0;
      for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = dataBase.nodeListBegin();
           nodeListItr != dataBase.nodeListEnd();
           ++nodeListItr, ++nodeListID) {
        ENSURE((int)(**nodeListItr).numNodes() == firstNewGhostNode[nodeListID]);
      }
    }
  }
  END_CONTRACT_SCOPE

}

}

