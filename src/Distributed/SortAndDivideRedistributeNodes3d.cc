//---------------------------------Spheral++----------------------------------//
// SortAndDivideRedistributeNodes -- (Re)domain decompose the nodes by using the
// information encoded in the SortAndDivideNeighbor algorithm.
//----------------------------------------------------------------------------//
#include "Distributed/SortAndDivideRedistributeNodes3d.hh"
#include "Distributed/TreeDistributedBoundary.hh"
#include "Utilities/DomainNode.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Geometry/EigenStruct.hh"
#include "Utilities/globalNodeIDs.hh"
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
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given node extent.
//------------------------------------------------------------------------------
SortAndDivideRedistributeNodes3d::
SortAndDivideRedistributeNodes3d(const double Hextent):
  SortAndDivideRedistributeNodes<Dim<3> >(Hextent) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SortAndDivideRedistributeNodes3d::
~SortAndDivideRedistributeNodes3d() {
}

//------------------------------------------------------------------------------
// The main method of this class.  Repartition the node distribution.
//------------------------------------------------------------------------------
void
SortAndDivideRedistributeNodes3d::
redistributeNodes(DataBase<Dim<3> >& dataBase,
                  vector<Boundary<Dim<3> >*> boundaries) {

  // Number of processors.
  const int procID = this->domainID();
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
  const EigenStruct<3> shapeTensor = this->shapeTensor(nodeDistribution);

  // Rotate the nodes into this frame.
  this->rotateIntoShapeTensorFrame(shapeTensor, nodeDistribution);

  // Get the number of domains we should have per constant work chunk in the primary direction.
  vector< vector<int> > domainsPerStep = domainsPerChunk(shapeTensor);

  // Output the initial load distribution statistics.
  const string initialLoadStats = this->gatherDomainDistributionStatistics(work);
  if (procID == 0) {
    cout << "SortAndDivideRedistributeNodes::redistributeNodes initial load balancing:" << endl
         << initialLoadStats << endl
         << "    Domain distribution shape tensor: " << shapeTensor.eigenValues << endl;
    for (int i = 0; i != Dimension::nDim; ++i) {
      cout << "    " << shapeTensor.eigenVectors.getColumn(i) << endl;
    }
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

    // Pop off the set of nodes for the y-z slab in this step.
    const vector<int>& numYZChunks = domainsPerStep[ix];
    CHECK(numYZChunks.size() > 0);
    int numYZDomains = 0;
    for (vector<int>::const_iterator itr = numYZChunks.begin();
         itr != numYZChunks.end();
         ++itr) numYZDomains += *itr;
    const double targetYZChunkWork = numYZDomains*targetWorkPerDomain;

    list<DomainNode<Dimension> > yzchunkNodes;
    if (ix == domainsPerStep.size() - 1) {
      yzchunkNodes = sortedNodeDistribution;
    } else {
      yzchunkNodes = this->popFrontNodes(sortedNodeDistribution,
                                         targetYZChunkWork,
                                         0);
    }

    // Re-sort this set of nodes by their y positions.
    this->sortByPositions(yzchunkNodes, 1);

    // Iterate over the number of y domains we're assigning for this slab of work.
    for (auto iy = 0u; iy != numYZChunks.size(); ++iy) {

      // Pop off a chunk of nodes we'll divvy up in the z direction.
      const int numZChunks = numYZChunks[iy];
      CHECK(numZChunks > 0);
      const double targetZChunkWork = numZChunks*targetWorkPerDomain;

      // If this is the last domain in this slice, make sure we got all the nodes.
      // Otherwise, just pop off the target domains' work.
      list<DomainNode<Dimension> > chunkNodes;
      if (iy == numYZChunks.size() - 1) {
        chunkNodes = yzchunkNodes;
      } else {
        chunkNodes = this->popFrontNodes(yzchunkNodes,
                                         targetZChunkWork,
                                         1);
      }

      // Re-sort this set of nodes by their z positions.
      this->sortByPositions(chunkNodes, 2);

      // Iterator over the number of z domains we'll be assigning.
      for (int iz = 0; iz != numZChunks; ++iz) {

        if (procID == 0) cout << "Assigning domain " << assignDomainID 
                              << " of " << numProcs << "...";

        // Peel off nodes from the front of the unassigned nodes, until the desired work
        // load for this domain is reached.  Note that in this step we use the z index as
        // the primary sorting comparison.
        list<DomainNode<Dimension> > thisDomainNodes;
        if (iz == numZChunks - 1) {
          thisDomainNodes = chunkNodes;
        } else {
          thisDomainNodes = this->popFrontNodes(chunkNodes,
                                                targetWorkPerDomain,
                                                2);
        }

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
        if (procID == 0) cout << "Done." << endl;

      }

//       // Assign any remaining nodes from this chunks work to the last domain we were working
//       // on.
//       CHECK(assignDomainID - 1 >= 0 && assignDomainID - 1 < numProcs);
//       for (list<DomainNode<Dimension> >::const_iterator itr = chunkNodes.begin();
//            itr != chunkNodes.end();
//            ++itr) {
//         newNodeDistribution.push_back(*itr);
//         newNodeDistribution.back().domainID = assignDomainID - 1;
//       }

    }

  }

//   // Assign any remaining nodes to the last domain.
//   for (list<DomainNode<Dimension> >::const_iterator itr = sortedNodeDistribution.begin();
//        itr != sortedNodeDistribution.end();
//        ++itr) {
//     newNodeDistribution.push_back(*itr);
//     newNodeDistribution.back().domainID = numProcs - 1;
//   }

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
  if (procID == 0) cout << "SortAndDivideRedistributeNodes::redistributeNodes final load balancing:" << endl
                        << finalLoadStats << endl << endl;
  MPI_Barrier(Communicator::communicator());

}

//------------------------------------------------------------------------------
// Compute the appropriate number of domains in each dimension for a given
// shape tensor EigenStruct.
//------------------------------------------------------------------------------
vector< vector<int> >
SortAndDivideRedistributeNodes3d::
domainsPerChunk(const Dim<3>::SymTensor::EigenStructType& shapeTensor) const {

  REQUIRE(fuzzyEqual(shapeTensor.eigenValues.magnitude(), 1.0));

  // The total number of domains we are going to assign.
  const int numProcs = this->numDomains();

  // Determine how many x-chunks we're going to use, and the corresponding
  // base number of domains we want in the yz slabs.
  const int xChunks = max(1, int(shapeTensor.eigenValues.x()/shapeTensor.eigenValues.sumElements() * numProcs));
  const int yzChunks = numProcs / xChunks;
  CHECK(xChunks*yzChunks > 0 && xChunks*yzChunks <= numProcs);

  // How many total remainder domains do we have to account for?
  const int totalRemainProcs = numProcs - xChunks*yzChunks;
  CHECK(totalRemainProcs >= 0);

  // Figure out how many remainder domains we have to account for in each
  // yz slab.
  const int baseRemainPerSlab = totalRemainProcs / xChunks;
  vector<int> remainProcs((std::size_t) xChunks, baseRemainPerSlab);
  remainProcs[0] += totalRemainProcs - baseRemainPerSlab * xChunks;
  BEGIN_CONTRACT_SCOPE
  {
    int checkCount = 0;
    CONTRACT_VAR(checkCount);
    for (const auto x: remainProcs) checkCount += x;
    CHECK(checkCount == totalRemainProcs);
  }
  END_CONTRACT_SCOPE

  // Prepare the result.
  vector< vector<int> > result(xChunks);

  // Iterate over each yz slab.
  for (int i = 0; i != xChunks; ++i) {
    CHECK(i < (int)result.size());
    CHECK(i < (int)remainProcs.size());

    // The total number of domains for this slab.
    const int numDomainsInSlab = yzChunks + remainProcs[i];

    // Determine how many y-chunks we're going to use, and the corresponding
    // base number of z-chunks.
    const int yChunks = max(1, int(shapeTensor.eigenValues.y()/(shapeTensor.eigenValues.y() + shapeTensor.eigenValues.z()) * numDomainsInSlab));
    const int zChunks = numDomainsInSlab / yChunks;
    CHECK(yChunks*zChunks > 0 && yChunks*zChunks <= numDomainsInSlab);

    // How many remainder domains do we have to account for?
    const int remainProcs = numDomainsInSlab - yChunks*zChunks;
    CHECK(remainProcs >= 0 && remainProcs < yChunks);

    // Build the numbers of domains in this slab.
    for (int j = 0; j != yChunks; ++j) {
      if (j < remainProcs) {
        result[i].push_back(zChunks + 1);
      } else {
        result[i].push_back(zChunks);
      }
    }
    CHECK((int)result[i].size() == yChunks);

    BEGIN_CONTRACT_SCOPE
    {
      int checkCount = 0;
      CONTRACT_VAR(checkCount);
      for (const auto x: result[i]) checkCount += x;
      CHECK(checkCount == numDomainsInSlab);
    }
    END_CONTRACT_SCOPE

  }

  // That's it.  Check our post-conditions and return the answer.
  ENSURE((int)result.size() == xChunks);
  BEGIN_CONTRACT_SCOPE
  {
    int checkCount = 0;
    CONTRACT_VAR(checkCount);
    for (int i = 0; i != xChunks; ++i) {
      for (const auto x: result[i]) checkCount += x;
    }
    ENSURE(checkCount == numProcs);
  }
  END_CONTRACT_SCOPE

  return result;
}

}
