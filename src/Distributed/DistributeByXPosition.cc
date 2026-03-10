//---------------------------------Spheral++----------------------------------//
// DistributeByXPosition -- Redistribute nodes by sorting their positions
// in x coordinate.  Really only useful in 1-D, as a test.
//
// Created by JMO, Mon Feb 10 16:12:47 PST 2003
//----------------------------------------------------------------------------//
#include <mpi.h>

#include "RedistributeNodes.hh"
#include "DistributeByXPosition.hh"
#include "Utilities/DomainNode.hh"
#include "DataBase/DataBase.hh"
#include "NodeList/NodeList.hh"
#include "Field/FieldList.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Communicator.hh"
#include "Utilities/DBC.hh"

#include <algorithm>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Local function to help sort a vector of DomainNode by x position.
//------------------------------------------------------------------------------
template<typename Dimension>
class xPositionLess {
public:
  bool operator()(const DomainNode<Dimension>& lhs,
                  const DomainNode<Dimension>& rhs) const {
  return lhs.position.x() < rhs.position.x();
  }
};

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
template<typename Dimension>
DistributeByXPosition<Dimension>::
DistributeByXPosition():
  RedistributeNodes<Dimension>() {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
DistributeByXPosition<Dimension>::
~DistributeByXPosition() {
}

//------------------------------------------------------------------------------
// The only real method of this class.  Sort the global positions of the nodes
// by x coordinate, and repartition the nodes across domains accordingly.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DistributeByXPosition<Dimension>::
redistributeNodes(DataBase<Dimension>& dataBase,
                  vector<Boundary<Dimension>*> /*boundaries*/) {

  // We're going to do this in a really dumb way (it's only a test for 1-D
  // anyway.)  Have all processors sort their own nodes positions by x,
  // then have process 0 collect the info from each process, and determine
  // the partitioning all by itself.

  // Assign temporary global IDs for each NodeList in the DataBase.
  const FieldList<Dimension, size_t> globalNodeIDs = Spheral::globalNodeIDs(dataBase);
  CHECK(globalNodeIDs.size() == dataBase.numNodeLists());

  // First get the local description of the domain distribution.
  vector<DomainNode<Dimension> > localDistribution = this->currentDomainDecomposition(dataBase, globalNodeIDs);

  // Sort the local localDistribution by x position.
  sort(localDistribution.begin(), localDistribution.end(), 
       xPositionLess<Dimension>());

  // Prepare a vector<DomainNode> to hold the global distribution.
  vector<DomainNode<Dimension> > globalDistribution = localDistribution;

  // Have each processor successively send it's results to processor 0, 
  // which accumulates the complete sorted list.
  const int procID = this->domainID();
  const int numProcs = this->numDomains();
  const int totalNumNodes = this->numGlobalNodes(dataBase);

  if (procID == 0) {

    // Process 0 receives from everyone else.
    for (int recvProc = 1; recvProc < numProcs; ++recvProc) {

      // Get the size of the buffer we're receiving.
      int sizeOfBuffer;
      {
        MPI_Status status;
        MPI_Recv(&sizeOfBuffer, 1, MPI_INT, recvProc, 11,
                 Communicator::communicator(), &status);
      }

      // Now get the encoded distribution from recvProc.
      CHECK(sizeOfBuffer >= 0);
      vector<char> recvBuffer(sizeOfBuffer);
      {
        MPI_Status status;
        MPI_Recv(&(*recvBuffer.begin()), sizeOfBuffer, MPI_CHAR, recvProc,
                 12, Communicator::communicator(), &status);
      }

      // Decode this message back into a set of DomainNode.
      vector<DomainNode<Dimension> > recvDistribution = this->unpackDomainNodes(recvBuffer);

      // Merge this set into the local info, in order of x position.
      vector<DomainNode<Dimension> > mergedResult(globalDistribution.size() +
                                                  recvDistribution.size());
      merge(globalDistribution.begin(), globalDistribution.end(),
            recvDistribution.begin(), recvDistribution.end(),
            mergedResult.begin(), xPositionLess<Dimension>());
      globalDistribution = mergedResult;
    }
    CHECK((int)globalDistribution.size() == totalNumNodes);

  } else {

    // All other processes just send their info to process 0.
    vector<char> encodedDistribution = this->packDomainNodes(localDistribution);
    int size = encodedDistribution.size();
    MPI_Send(&size, 1, MPI_INT, 0, 11, Communicator::communicator());
    MPI_Send(&(*encodedDistribution.begin()), size, MPI_CHAR, 0, 12, Communicator::communicator());
  }

  // OK, now process 0 has the complete set of sorted positions.  It will
  // now just divide this up between the other processors, and broadcast
  // the results to everyone else.
  if (procID == 0) {

    // Figure out how many nodes we want per domain.
    const int nominalNum = totalNumNodes/numProcs;
    const int remainder = totalNumNodes - nominalNum*numProcs;
    CHECK(remainder < numProcs);
    CHECK(nominalNum*numProcs + remainder == totalNumNodes);
    vector<int> numNodesPerDomain(numProcs);
    for (int i = 0; i < numProcs; ++i) {
      if (i < remainder) {
        numNodesPerDomain[i] = nominalNum + 1;
      } else {
        numNodesPerDomain[i] = nominalNum;
      }
    }

    // Now build the new domain distribution.
    int dID = 0;
    int lastStep = 0;
    for (auto i = 0u; i < globalDistribution.size(); ++i) {
      if ((int)i == lastStep + numNodesPerDomain[dID]) {
        lastStep = lastStep + numNodesPerDomain[dID];
        dID++;
      }
      CHECK(dID < numProcs);
      globalDistribution[i].domainID = dID;
    }
  }

  // Have process 0 broadcast the results to everyone.
  vector<char> encodedDistribution = this->packDomainNodes(globalDistribution);
  int size = encodedDistribution.size();
  MPI_Bcast(&size, 1, MPI_INT, 0, Communicator::communicator());
  if (procID > 0) encodedDistribution.resize(size);
  CHECK((int)encodedDistribution.size() == size);
  MPI_Bcast(&(*encodedDistribution.begin()), size, MPI_CHAR, 0, Communicator::communicator());
  if (procID > 0) globalDistribution = this->unpackDomainNodes(encodedDistribution);
  CHECK((int)globalDistribution.size() == totalNumNodes);

  // Now, reconstruct the local domain distribution with the new domain assignments.
  typename vector<DomainNode<Dimension> >::const_iterator globalItr = globalDistribution.begin();
  for (typename vector<DomainNode<Dimension> >::iterator localItr = localDistribution.begin();
       localItr < localDistribution.end();
       ++localItr) {
    while(globalItr < globalDistribution.end() &&
          globalItr->globalNodeID != localItr->globalNodeID) ++globalItr;
    CHECK(globalItr != globalDistribution.end());
    localItr->domainID = globalItr->domainID;
  }
  CHECK(globalItr != globalDistribution.end());

  // OK, the localDistribution now holds the desired redistribution of the nodes.
  // Go ahead and redistribute them.
  this->enforceDomainDecomposition(localDistribution, dataBase);
}

}

