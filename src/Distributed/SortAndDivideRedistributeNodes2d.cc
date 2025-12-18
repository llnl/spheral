//---------------------------------Spheral++----------------------------------//
// SortAndDivideRedistributeNodes -- (Re)domain decompose the nodes by using the
// information encoded in the SortAndDivideNeighbor algorithm.
//
// Created by JMO, Wed Nov 24 10:51:32 2004
//----------------------------------------------------------------------------//
#include "Distributed/SortAndDivideRedistributeNodes2d.hh"
#include "Distributed/TreeDistributedBoundary.hh"
#include "Utilities/DomainNode.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Geometry/EigenStruct.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/SpheralMessage.hh"
#include "Communicator.hh"

#include "Utilities/DBC.hh"

#include <algorithm>
#include <vector>
#include <map>
#include <fstream>
#include <cstdlib>
using std::vector;
using std::list;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given node extent.
//------------------------------------------------------------------------------
SortAndDivideRedistributeNodes2d::
SortAndDivideRedistributeNodes2d(const double Hextent):
  SortAndDivideRedistributeNodes<Dim<2> >(Hextent) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SortAndDivideRedistributeNodes2d::
~SortAndDivideRedistributeNodes2d() {
}

//------------------------------------------------------------------------------
// The main method of this class.  Repartition the node distribution.
//------------------------------------------------------------------------------
void
SortAndDivideRedistributeNodes2d::
redistributeNodes(DataBase<Dim<2> >& dataBase,
                  vector<Boundary<Dim<2> >*> boundaries) {

  // Number of processors.
  const int numProcs = this->numDomains();

  // Go over each NodeList, and clear out any ghost nodes.
  for (DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    (*nodeListItr)->numGhostNodes(0);
    (*nodeListItr)->neighbor().updateNodes();
  }

  // If the user did not specify any boundaries, then create a Distributed
  // boundary for local use.
  // Note that if boundary conditions were passed in, we assume that the Distributed
  // boundary is already in there.
  TreeDistributedBoundary<Dimension>& bound = TreeDistributedBoundary<Dimension>::instance();
  if (boundaries.size() == 0) boundaries.push_back(&bound);

  // Use the boundaries to create the ghost nodes.  We need this to compute the
  // work per node (defined as the number of neighbors) correctly.
  for (vector<Boundary<Dimension>*>::iterator boundItr = boundaries.begin();
       boundItr != boundaries.end();
       ++boundItr) {
    (*boundItr)->setAllGhostNodes(dataBase);
    (*boundItr)->finalizeGhostBoundary();
    for (DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
         nodeListItr != dataBase.nodeListEnd();
         ++nodeListItr) (*nodeListItr)->neighbor().updateNodes();
  }

  // Get the work per node.
  const FieldList<Dimension, Scalar> work = this->workPerNode(dataBase, Hextent());

  // Once again clear out any ghost nodes.
  for (DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    (*nodeListItr)->numGhostNodes(0);
    (*nodeListItr)->neighbor().updateNodes();
  }

  // Build the set of global node IDs.
  const FieldList<Dimension, size_t> globalIDs = globalNodeIDs(dataBase);

  // Get the local description of the domain distribution.
  vector<DomainNode<Dimension> > nodeDistribution = this->currentDomainDecomposition(dataBase, globalIDs, work);

  // Clear the domain assignments of all nodes.
  for (vector<DomainNode<Dimension> >::iterator itr = nodeDistribution.begin();
       itr != nodeDistribution.end();
       ++itr) itr->domainID = -1;

  // Compute the shape tensor describing the total node distribution in space.
  const EigenStruct<2> shapeTensor = this->shapeTensor(nodeDistribution);

  // Rotate the nodes into this frame.
  this->rotateIntoShapeTensorFrame(shapeTensor, nodeDistribution);

  // Get the number of domains we should have per constant work chunk in the primary direction.
  vector<int> domainsPerStep = domainsPerChunk(shapeTensor);

  // Output the initial load distribution statistics.
  const string initialLoadStats = this->gatherDomainDistributionStatistics(work);
  if (Process::getRank() == 0) {
    SpheralMessage("SortAndDivideRedistributeNodes::redistributeNodes initial load balancing:\n" << 
                   initialLoadStats << "\n" <<
                   "    Domain distribution shape tensor: " << shapeTensor.eigenValues << "\n" <<
                                                 "    Number of domains per work chunk: ");
    for (auto& x: domainsPerStep) SpheralMessage("                                      " << x);
  }
    
  // Compute the total work, and the target work per processor.
  double localWork = 0.0;
  for (vector<DomainNode<Dimension> >::const_iterator itr = nodeDistribution.begin();
       itr != nodeDistribution.end();
       ++itr) localWork += itr->work;
  CHECK(numProcs > 0);
  double globalWork;
  MPI_Allreduce(&localWork, &globalWork, 1, MPI_DOUBLE, MPI_SUM, Communicator::communicator());
  const double targetWorkPerDomain = globalWork / numProcs;
  CHECK(distinctlyGreaterThan(targetWorkPerDomain, 0.0));

  // Copy the domain nodes to a list.
  list<DomainNode<Dimension> > sortedNodeDistribution;
  for (vector<DomainNode<Dimension> >::const_iterator itr = nodeDistribution.begin();
       itr != nodeDistribution.end();
       ++itr) sortedNodeDistribution.push_back(*itr);

  // Sort the nodes by x position.
  this->sortByPositions(sortedNodeDistribution, 0);

  // The new node distribution we're going to build.
  vector<DomainNode<Dimension> > newNodeDistribution;
  newNodeDistribution.reserve(nodeDistribution.size());

  // Iterate over the x domain indicies.
  int assignDomainID = 0;
  for (auto ix = 0u; ix != domainsPerStep.size(); ++ix) {

    // Pop off the set of nodes we'll divvy up in the y direction in this step.
    const int numYChunks = domainsPerStep[ix];
    CHECK(numYChunks > 0);
    const double targetChunkWork = numYChunks*targetWorkPerDomain;
    list<DomainNode<Dimension> > chunkNodes = this->popFrontNodes(sortedNodeDistribution,
                                                                  targetChunkWork,
                                                                  0);

    // Re-sort this set of nodes by their y positions.
    this->sortByPositions(chunkNodes, 1);

    // Iterate over the number of y domains we're assigning for this x chunk of work.
    for (int iy = 0; iy != numYChunks; ++iy) {

      // Peel off nodes from the front of the unassigned nodes, until the desired work
      // load for this domain is reached.  Note that in this step we use the y index as
      // the primary sorting comparison.
      const list<DomainNode<Dimension> > thisDomainNodes = this->popFrontNodes(chunkNodes,
                                                                               targetWorkPerDomain,
                                                                               1);

      // Assign the nodes to their new domain.
      CHECK(assignDomainID >= 0 && assignDomainID < numProcs);
      for (list<DomainNode<Dimension> >::const_iterator itr = thisDomainNodes.begin();
           itr != thisDomainNodes.end();
           ++itr) {
        newNodeDistribution.push_back(*itr);
        newNodeDistribution.back().domainID = assignDomainID;
      }

      // Increment the domain we're assigning to.
      ++assignDomainID;

    }

    // Assign any remaining nodes from this chunks work to the last domain we were working
    // on.
    CHECK(assignDomainID - 1 >= 0 && assignDomainID - 1 < numProcs);
    for (list<DomainNode<Dimension> >::const_iterator itr = chunkNodes.begin();
         itr != chunkNodes.end();
         ++itr) {
      newNodeDistribution.push_back(*itr);
      newNodeDistribution.back().domainID = assignDomainID - 1;
    }

  }

  // Assign any remaining nodes to the last domain.
  for (list<DomainNode<Dimension> >::const_iterator itr = sortedNodeDistribution.begin();
       itr != sortedNodeDistribution.end();
       ++itr) {
    newNodeDistribution.push_back(*itr);
    newNodeDistribution.back().domainID = numProcs - 1;
  }

  // Check that all nodes really really have been assigned to a domain.
  BEGIN_CONTRACT_SCOPE
  CHECK(newNodeDistribution.size() == nodeDistribution.size());
  for (vector<DomainNode<Dimension> >::iterator itr = newNodeDistribution.begin();
       itr != newNodeDistribution.end();
       ++itr) CHECK(itr->domainID >= 0 && itr->domainID < numProcs);
  END_CONTRACT_SCOPE

  // The nodeDistribution now holds the desired redistribution of the nodes.
  // Go ahead and redistribute them.
  enforceDomainDecomposition(newNodeDistribution, dataBase);

  // Output the final load distribution statistics.
  const string finalLoadStats = this->gatherDomainDistributionStatistics(work);
  SpheralMessage("SortAndDivideRedistributeNodes::redistributeNodes final load balancing:\n" << finalLoadStats << "\n\n");
  MPI_Barrier(Communicator::communicator());
}

//------------------------------------------------------------------------------
// Compute the appropriate number of domains in each dimension for a given
// shape tensor EigenStruct.
//------------------------------------------------------------------------------
std::vector<int> 
SortAndDivideRedistributeNodes2d::
domainsPerChunk(const Dim<2>::SymTensor::EigenStructType& shapeTensor) const {

  REQUIRE(fuzzyEqual(shapeTensor.eigenValues.magnitude(), 1.0));

  // The total number of domains we are going to assign.
  const int numProcs = this->numDomains();

  // Determine how many x-chunks we're going to use, and the corresponding
  // base number of y-chunks.
  const int xChunks = max(1, int(shapeTensor.eigenValues.x()/shapeTensor.eigenValues.sumElements() * numProcs));
  const int yChunks = numProcs / xChunks;
  CHECK(xChunks*yChunks > 0 && xChunks*yChunks <= numProcs);

  // How many remainder domains do we have to account for?
  const int remainProcs = numProcs - xChunks*yChunks;
  CHECK(remainProcs >= 0 && remainProcs < xChunks);

  // Now we can build the result.
  vector<int> result(xChunks);
  for (int i = 0; i != xChunks; ++i) {
    if (i < remainProcs) {
      result[i] = yChunks + 1;
    } else {
      result[i] = yChunks;
    }
  }

  // That's it.  Check our post-conditions and return the answer.
  BEGIN_CONTRACT_SCOPE
  {
    int checkCount = 0;
    CONTRACT_VAR(checkCount);
    for (const auto x: result) checkCount += x;
    ENSURE(checkCount == numProcs);
  }
  END_CONTRACT_SCOPE

  return result;
}

}
