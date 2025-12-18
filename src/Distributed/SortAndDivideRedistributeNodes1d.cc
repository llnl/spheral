//---------------------------------Spheral++----------------------------------//
// SortAndDivideRedistributeNodes -- (Re)domain decompose the nodes by using the
// information encoded in the SortAndDivideNeighbor algorithm.
//
// Created by JMO, Wed Nov 24 10:51:32 2004
//----------------------------------------------------------------------------//
#include "Distributed/SortAndDivideRedistributeNodes1d.hh"
#include "Distributed/TreeDistributedBoundary.hh"
#include "Utilities/DomainNode.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
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
SortAndDivideRedistributeNodes1d::
SortAndDivideRedistributeNodes1d(const double Hextent):
  SortAndDivideRedistributeNodes<Dim<1> >(Hextent) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SortAndDivideRedistributeNodes1d::
~SortAndDivideRedistributeNodes1d() {
}

//------------------------------------------------------------------------------
// The main method of this class.  Repartition the node distribution.
//------------------------------------------------------------------------------
void
SortAndDivideRedistributeNodes1d::
redistributeNodes(DataBase<Dim<1> >& dataBase,
                  vector<Boundary<Dim<1> >*> boundaries) {

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

  // Output the initial load distribution statistics.
  const string initialLoadStats = this->gatherDomainDistributionStatistics(work);
  SpheralMessage("SortAndDivideRedistributeNodes::redistributeNodes initial load balancing:\n" << initialLoadStats << "\n\n");

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

  // Iterate over each domain we're going to assign to.
  for (int assignDomainID = 0; assignDomainID != numProcs - 1; ++assignDomainID) {

    // Peel off nodes from the front of the unassigned nodes, until the desired work
    // load is reached.
    const list<DomainNode<Dimension> > thisDomainNodes = this->popFrontNodes(sortedNodeDistribution,
                                                                             targetWorkPerDomain,
                                                                             0);

    // Assign the nodes to their new domain.
    for (list<DomainNode<Dimension> >::const_iterator itr = thisDomainNodes.begin();
         itr != thisDomainNodes.end();
         ++itr) {
      newNodeDistribution.push_back(*itr);
      newNodeDistribution.back().domainID = assignDomainID;
    }

  }

  // Now there is just one domain left.  Assign the remaining nodes to it.
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
  CHECK(validDomainDecomposition(newNodeDistribution, dataBase));
  enforceDomainDecomposition(newNodeDistribution, dataBase);

  // Output the final load distribution statistics.
  const string finalLoadStats = this->gatherDomainDistributionStatistics(work);
  SpheralMessage("SortAndDivideRedistributeNodes::redistributeNodes final load balancing:\n" << finalLoadStats << "\n\n");
  MPI_Barrier(Communicator::communicator());

}

}
