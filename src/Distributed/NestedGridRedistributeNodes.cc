//---------------------------------Spheral++----------------------------------//
// NestedGridRedistributeNodes -- (Re)domain decompose the nodes by using the
// information encoded in the NestedGridNeighbor algorithm.
//
// Created by JMO, Wed Nov 24 10:51:32 2004
//----------------------------------------------------------------------------//
#include "NestedGridRedistributeNodes.hh"
#include "Utilities/DomainNode.hh"
#include "NestedGridDistributedBoundary.hh"
#include "DistributedBoundary.hh"
#include "NestedGridUtilities.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/NestedGridNeighbor.hh"
#include "Neighbor/GridCellIndex.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/SpheralMessage.hh"
#include "Communicator.hh"

#include "Utilities/DBC.hh"

#include <algorithm>
#include <vector>
#include <map>

#include <fstream>
#include <cstdlib>


namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given node extent.
//------------------------------------------------------------------------------
template<typename Dimension>
NestedGridRedistributeNodes<Dimension>::
NestedGridRedistributeNodes(const double Hextent):
  mHextent(Hextent) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
NestedGridRedistributeNodes<Dimension>::
~NestedGridRedistributeNodes() {
}

//------------------------------------------------------------------------------
// The main method of this class.  Repartition the node distribution.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridRedistributeNodes<Dimension>::
redistributeNodes(DataBase<Dimension>& dataBase,
                  vector<Boundary<Dimension>*> boundaries) {

  // Number of processors.
  const int procID = this->domainID();
  const int numProcs = this->numDomains();

  // Go over each NodeList, and clear out any ghost nodes.
  for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    (*nodeListItr)->numGhostNodes(0);
    (*nodeListItr)->neighbor().updateNodes();
  }

  // Build the set of global node IDs.
  const FieldList<Dimension, int> globalIDs = globalNodeIDs(dataBase);

  // Get the local description of the domain distribution.
  vector<DomainNode<Dimension> > nodeDistribution = this->currentDomainDecomposition(dataBase, globalIDs);
  const int numLocalNodes = nodeDistribution.size();

  // If the user did not specify any boundaries, then create a Distributed
  // boundary for local use.
  // Note that if boundary conditions were passed in, we assume that the Distributed
  // boundary is already in there.
  NestedGridDistributedBoundary<Dimension>& bound = 
    NestedGridDistributedBoundary<Dimension>::instance();
  if (boundaries.size() == 0) boundaries.push_back(&bound);

  // Use the boundaries to create the ghost nodes.  We need this to compute the
  // work per node (defined as the number of neighbors) correctly.
  for (typename vector<Boundary<Dimension>*>::iterator boundItr = boundaries.begin();
       boundItr != boundaries.end();
       ++boundItr) {
    (*boundItr)->setAllGhostNodes(dataBase);
    (*boundItr)->finalizeGhostBoundary();
    for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
         nodeListItr != dataBase.nodeListEnd();
         ++nodeListItr) (*nodeListItr)->neighbor().updateNodes();
  }

  // Get the work per node.
  const FieldList<Dimension, Scalar> work = this->workPerNode(dataBase,
                                                              mHextent);

  // Once again clear out any ghost nodes.
  for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    (*nodeListItr)->numGhostNodes(0);
    (*nodeListItr)->neighbor().updateNodes();
  }

  // Output the initial load distribution statistics.
  const string initialLoadStats = this->gatherDomainDistributionStatistics(work);
  SpheralMessage("NestedGridRedistributeNodes::redistributeNodes initial load balancing:\n" << initialLoadStats << "\n\n");

  // Compute the total work, and the target work per processor.
  double localWork = 0.0;
  for (int i = 0; i != work.size(); ++i) {
    localWork += work[i]->sumElements();
  }
  CHECK(numProcs > 0);
  double globalWork;
  MPI_Allreduce(&localWork, &globalWork, 1, MPI_DOUBLE, MPI_SUM, Communicator::communicator());
  const double targetWorkPerDomain = globalWork / numProcs;
  CHECK(distinctlyGreaterThan(targetWorkPerDomain, 0.0));

  // Get the global counts of nodes in each grid cell and each grid level.
  GridCellPopulationType gridCellPopulations;
  vector<int> gridLevelPopulations;
  computeGridCellPopulations(dataBase, gridCellPopulations, gridLevelPopulations);
  CHECK(gridCellPopulations.size() == gridLevelPopulations.size());

  // Clear the domain assignments of all nodes.
  for (typename vector<DomainNode<Dimension> >::iterator itr = nodeDistribution.begin();
       itr != nodeDistribution.end();
       ++itr) itr->domainID = -1;

  // Iterate until all nodes have been assigned to a domain.
  int currentDomain = 0;
  Scalar currentDomainWork = 0.0;
  bool currentDomainFull = false;
  while (countRemainingNodes(gridLevelPopulations) > 0) {
    CHECK(countRemainingNodes(gridLevelPopulations) == countRemainingNodes(gridCellPopulations));
    CHECK(currentDomain >= 0 && currentDomain < numProcs);
//     cerr << "Remaining nodes: " << countRemainingNodes(gridLevelPopulations) << endl;

    // Find the grid level with the greatest number of remaining nodes, and the minimum
    // remaining grid cell on that grid level.
    GridCellIndex<Dimension> masterGridCell;
    int masterGridLevel;
    findMaster(gridCellPopulations, gridLevelPopulations, masterGridCell, masterGridLevel);

    // Fill the current domain.
    int gridCellStep = 0;
    while (!currentDomainFull && countRemainingNodes(gridLevelPopulations) > 0) {

      // Check to see if there are still nodes left on the current master grid level.
      if (gridLevelPopulations[masterGridLevel] == 0) {

        findMaster(gridCellPopulations, gridLevelPopulations, masterGridCell, masterGridLevel);
        gridCellStep = 0;

//         // The old master grid level is empty!  We need to pick a new master grid
//         // level, trying to keep to the same spatial location.
//         int newGridLevel;
//         GridCellIndex<Dimension> minGridCell;
//         findMaster(gridCellPopulations, gridLevelPopulations, minGridCell, newGridLevel);
//         const int dgl = abs(newGridLevel - masterGridLevel);
//         const int scaleFactor = int(pow(2.0, dgl) + 0.1);
//         CHECK(scaleFactor > 0);
//         if (newGridLevel > masterGridLevel) {
//           masterGridCell = masterGridCell*scaleFactor;
//         } else {
//           masterGridCell = masterGridCell/scaleFactor;
//         }
//         gridCellStep = 0;

      }

      // Create the potential set of master grid cells.
      const vector<GridCellIndex<Dimension> > masterGridCells = gridCellRind(masterGridCell,
                                                                             masterGridLevel,
                                                                             gridCellStep,
                                                                             gridCellPopulations);
      ++gridCellStep;

      // Iterate over the set of potential master grid cells, until we either exhaust that set or
      // have filled the current domain.
      typename vector<GridCellIndex<Dimension> >::const_iterator gcItr = masterGridCells.begin();
//       cerr << (gcItr != masterGridCells.end()) << " " << (!currentDomainFull) << " " << countRemainingNodes(gridLevelPopulations) << endl;
      while (gcItr != masterGridCells.end() && 
             !currentDomainFull && 
             countRemainingNodes(gridLevelPopulations) > 0) {
//         cerr <<   "Checking grid cell/grid cell step: " << *gcItr << " " << gridCellStep << endl;

        // Set the master information for all NodeLists (and therefore coarse neighbors as well).
        vector<vector<int>> localMasterLists(dataBase.numNodeLists()), localCoarseNeighbors(dataBase.numNodeLists());
        setMasterNodeLists(dataBase, *gcItr, masterGridLevel, localMasterLists, localCoarseNeighbors);

        // Gather the available nodes from the coarse neighbor set to process 0.
        vector<int> coarseNodeIndices;
        vector<Scalar> coarseNodeWork;
        gatherAvailableCoarseNodes(dataBase, nodeDistribution, work, localCoarseNeighbors, coarseNodeIndices, coarseNodeWork);

        // Assign nodes from this set to the current domain, until it's work quota is
        // filled.
        currentDomainFull = assignNodesToDomain(dataBase,
                                                coarseNodeIndices,
                                                coarseNodeWork,
                                                currentDomain,
                                                targetWorkPerDomain,
                                                currentDomainWork,
                                                nodeDistribution,
                                                gridCellPopulations,
                                                gridLevelPopulations);

        // Increment the current grid cell.
        ++gcItr;
      }
    }

    // Go on to a new domain.
    ++currentDomain;
    currentDomainWork = 0.0;
    currentDomainFull = false;
  }
  CHECK(currentDomain == numProcs);

  // Check that all nodes really really have been assigned to a domain.
  BEGIN_CONTRACT_SCOPE
  for (typename vector<DomainNode<Dimension> >::iterator itr = nodeDistribution.begin();
       itr != nodeDistribution.end();
       ++itr) CHECK(itr->domainID >= 0 && itr->domainID < numProcs);
  END_CONTRACT_SCOPE

  // The nodeDistribution now holds the desired redistribution of the nodes.
  // Go ahead and redistribute them.
  CHECK(this->validDomainDecomposition(nodeDistribution, dataBase));
  this->enforceDomainDecomposition(nodeDistribution, dataBase);

  // Output the final load distribution statistics.
  const string finalLoadStats = this->gatherDomainDistributionStatistics(work);
  SpheralMessage("NestedGridRedistributeNodes::redistributeNodes final load balancing:\n" << finalLoadStats << "\n\n");
  MPI_Barrier(Communicator::communicator());

}

//------------------------------------------------------------------------------
// Compute the number of nodes in each grid cell, and the total number on
// each grid level.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridRedistributeNodes<Dimension>::
computeGridCellPopulations(const DataBase<Dimension>& dataBase,
                           vector< map<GridCellIndex<Dimension>, int> >& gridCellPopulations,
                           vector<int>& gridLevelPopulations) const {

  // Make sure the return variables are sized properly and zeroed.
  const size_t numOccupiedGridLevels = maxNumGridLevels(dataBase, Communicator::communicator());
  gridCellPopulations = GridCellPopulationType(numOccupiedGridLevels);
  gridLevelPopulations = vector<int>(numOccupiedGridLevels, 0);

  // Iterate over the NodeLists in the DataBase.
  for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {

    // Get the set of occupied grid cells for this NodeList.
    const NestedGridNeighbor<Dimension>& neighbor = getNestedGridNeighbor(*nodeListItr);
    const vector< vector<GridCellIndex<Dimension> > >& occupiedGridCells = neighbor.occupiedGridCells();

    // Iterate over the grid levels.
    for (int gridLevel = 0; gridLevel != numOccupiedGridLevels; ++gridLevel) {
      CHECK(gridLevel < occupiedGridCells.size());
      CHECK(gridLevel < gridCellPopulations.size());
      CHECK(gridLevel < gridLevelPopulations.size());

      // Iterate over the occupied grid cells on this grid level for this NodeList.
      for (typename vector<GridCellIndex<Dimension> >::const_iterator gcItr = occupiedGridCells[gridLevel].begin();
           gcItr != occupiedGridCells[gridLevel].end();
           ++gcItr) {
        const int numNodesInGridCell = neighbor.nodesInCell(*gcItr, gridLevel).size();

        // Increment the counts of nodes in this grid cell and on this grid level.
        if (gridCellPopulations[gridLevel].find(*gcItr) == gridCellPopulations[gridLevel].end()) gridCellPopulations[gridLevel][*gcItr] = 0;
        gridCellPopulations[gridLevel][*gcItr] += numNodesInGridCell;
        gridLevelPopulations[gridLevel] += numNodesInGridCell;

      }
    }
  }

  BEGIN_CONTRACT_SCOPE
  // Check that the local grid cell and grid level counts are consistent before
  // we begin to reduce global info.
  for (int gridLevel = 0; gridLevel != numOccupiedGridLevels; ++gridLevel) {
    int checkCount = 0;
    for (typename map<GridCellIndex<Dimension>, int>::const_iterator itr = gridCellPopulations[gridLevel].begin();
         itr != gridCellPopulations[gridLevel].end();
         ++itr) checkCount += itr->second;
    CHECK(checkCount == gridLevelPopulations[gridLevel]);
  }
  END_CONTRACT_SCOPE

  // Every domain has individually counted up the number of nodes, so sum the
  // counts over all processes to get the global counts.
  MPI_Barrier(Communicator::communicator());
  const int numProcs = this->numDomains();
  const int procID = this->domainID();
  for (int gridLevel = 0; gridLevel != numOccupiedGridLevels; ++gridLevel) {

    // Reduce the total numbers of nodes on the grid levels.
    int globalNumNodes;
    MPI_Allreduce(&(gridLevelPopulations[gridLevel]), &globalNumNodes, 1, MPI_INT, MPI_SUM, Communicator::communicator());
    gridLevelPopulations[gridLevel] = globalNumNodes;

    // The complication with grid cells is that each domain may have different
    // sets of occupied grid cells.  Therefore we have to have each process broadcast
    // its full information to all others.

    // Make a copy of the local grid cell population on this grid level.
    const map<GridCellIndex<Dimension>, int> gridCellPopulationsCopy(gridCellPopulations[gridLevel]);

    for (int sendProc = 0; sendProc != numProcs; ++sendProc) {

      if (procID == sendProc) {

        // Pack up our grid cell info.
        vector<int> buffer;
        for (typename map<GridCellIndex<Dimension>, int>::const_iterator itr = gridCellPopulationsCopy.begin();
             itr != gridCellPopulationsCopy.end();
             ++itr) {
          for (int i = 0; i != Dimension::nDim; ++i) buffer.push_back(itr->first(i));
          buffer.push_back(itr->second);
        }
        int bufferSize = buffer.size();
        CHECK(bufferSize == gridCellPopulationsCopy.size() * (Dimension::nDim + 1));

        // Send to all other processors.
        for (int proc = 0; proc != numProcs; ++proc) {
          if (proc != procID) {
            MPI_Send(&bufferSize, 1, MPI_INT, proc, 110, Communicator::communicator());
            MPI_Send(&(*buffer.begin()), bufferSize, MPI_INT, proc, 111, Communicator::communicator());
          }
        }

      } else {

        // Get the current send procs encoded info.
        int bufferSize;
        MPI_Status status1, status2;
        MPI_Recv(&bufferSize, 1, MPI_INT, sendProc, 110, Communicator::communicator(), &status1);
        vector<int> buffer(bufferSize);
        MPI_Recv(&(*buffer.begin()), bufferSize, MPI_INT, sendProc, 111, Communicator::communicator(), &status2);

        // Unpack the grid cell info and add it to our current set.
        typename vector<int>::const_iterator bufferItr = buffer.begin();
        while (bufferItr < buffer.end()) {
          GridCellIndex<Dimension> gc;
          for (int i = 0; i != Dimension::nDim; ++i) {
            CHECK(bufferItr < buffer.end());
            gc(i) = *bufferItr;
            ++bufferItr;
          }
          CHECK(bufferItr < buffer.end());
          const int numNodes = *bufferItr;
          ++bufferItr;

          if (gridCellPopulations[gridLevel].find(gc) == gridCellPopulations[gridLevel].end()) gridCellPopulations[gridLevel][gc] = 0;
          gridCellPopulations[gridLevel][gc] += numNodes;
        }
      }

    }
  }

  // Post-conditions.
  ENSURE(countRemainingNodes(gridLevelPopulations) == countRemainingNodes(gridCellPopulations));
}

//------------------------------------------------------------------------------
// Find the grid level that has the most nodes, and the minimum grid cell on 
// that grid level.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridRedistributeNodes<Dimension>::
findMaster(const typename NestedGridRedistributeNodes<Dimension>::GridCellPopulationType& gridCellPopulations,
           const vector<int>& gridLevelPopulations,
           GridCellIndex<Dimension>& gridCell,
           int& gridLevel) const {

  // Preconditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(gridCellPopulations.size() == gridLevelPopulations.size());
  for (int i = 0; i != gridCellPopulations.size(); ++i) {
    int check = 0;
    for (typename map<GridCellIndex<Dimension>, int>::const_iterator itr = gridCellPopulations[i].begin();
         itr != gridCellPopulations[i].end();
         ++itr) check += itr->second;
    REQUIRE(check == gridLevelPopulations[i]);
  }
  END_CONTRACT_SCOPE

  // First find the grid level with the most nodes.
  const int numGridLevels = gridLevelPopulations.size();
  gridLevel = 0;
  int maxGridLevelCount = 0;
  for (int i = 0; i != numGridLevels; ++i) {
    if (gridLevelPopulations[i] > maxGridLevelCount) {
      gridLevel = i;
      maxGridLevelCount = gridLevelPopulations[i];
    }
  }
  CHECK(gridLevel >= 0 && gridLevel < numGridLevels);

  // Now find the minimum (as defined by grid cell comparison) grid cell on this
  // grid level.
  CHECK(gridCellPopulations[gridLevel].size() > 0);
  gridCell = gridCellPopulations[gridLevel].begin()->first;
  int maxGridCellCount = 0;
  for (typename map<GridCellIndex<Dimension>, int>::const_iterator itr = gridCellPopulations[gridLevel].begin();
       itr != gridCellPopulations[gridLevel].end();
       ++itr) {
    if (itr->second > maxGridCellCount) {
      gridCell = itr->first;
      maxGridCellCount = itr->second;
    }
  }

  // That's it.  Both the grid cell and grid level should be set.
  ENSURE(gridLevel >= 0 && gridLevel < numGridLevels);
}

//------------------------------------------------------------------------------
// Count the number of nodes in the given grid cell per grid level population.
//------------------------------------------------------------------------------
template<typename Dimension>
int
NestedGridRedistributeNodes<Dimension>::
countRemainingNodes(const vector< map<GridCellIndex<Dimension>, int> >& gridCellPopulations) const {

  // Prepare the result.
  int result = 0;

  // Iterate over each grid level.
  for (typename GridCellPopulationType::const_iterator glItr = gridCellPopulations.begin();
       glItr != gridCellPopulations.end();
       ++glItr) {

    // Iterate over the grid cell entries for this grid level.
    for (typename map<GridCellIndex<Dimension>, int>::const_iterator gcItr = glItr->begin();
         gcItr != glItr->end();
         ++gcItr) {

      // Increment by the number of nodes in this grid cell.
      result += gcItr->second;

    }
  }

  return result;
}

//------------------------------------------------------------------------------
// Count the number of nodes in the given grid level population.
//------------------------------------------------------------------------------
template<typename Dimension>
int
NestedGridRedistributeNodes<Dimension>::
countRemainingNodes(const vector<int>& gridLevelPopulations) const {

  // Just sum up the counts on each grid level.
  int result = 0;
  for (typename vector<int>::const_iterator glItr = gridLevelPopulations.begin();
       glItr != gridLevelPopulations.end();
       ++glItr) result += *glItr;
  return result;
}

//------------------------------------------------------------------------------
// Select the set of populated grid cells a given spacing out from a specified
// center grid cell.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<GridCellIndex<Dimension> >
NestedGridRedistributeNodes<Dimension>::
gridCellRind(const GridCellIndex<Dimension>& gridCell,
             const int gridLevel,
             const int gridCellStep,
             const typename NestedGridRedistributeNodes::GridCellPopulationType& gridCellPopulations) const {

  // Pre-conditions.
  REQUIRE(gridLevel >= 0 && gridLevel < gridCellPopulations.size());
  REQUIRE(gridCellStep >= 0);

  // Prepare the result.
  vector<GridCellIndex<Dimension> > result;

  // Estimate how many grid cells we would have to walk if we just chunked over 
  // every one at the specified step.  This estimate is high by some duplicates, but 
  // is good enough for our purposes.
  const int numGridCellsAtStep = 2*Dimension::nDim * int(Dimension::pownu1(double(1 + 2*gridCellStep)) + 0.1);
  CHECK(numGridCellsAtStep > 0);

  // Now is it more efficient to walk every grid cell at the specified distance, or 
  // every occupied grid cell on the grid level?
  if (numGridCellsAtStep < gridCellPopulations[gridLevel].size()) {

    // It's cheaper to explictly walk the rind, so do it.
    const set<GridCellIndex<Dimension> > potentialGridCells = computeGridCellRind(gridCell, gridCellStep);
    CHECK(potentialGridCells.size() <= numGridCellsAtStep);
    CHECK(potentialGridCells.size() < gridCellPopulations[gridLevel].size());
    for (typename set<GridCellIndex<Dimension> >::const_iterator gcItr = potentialGridCells.begin();
         gcItr != potentialGridCells.end();
         ++gcItr) {
      if (gridCellPopulations[gridLevel].find(*gcItr) != gridCellPopulations[gridLevel].end()) result.push_back(*gcItr);
    }

  } else {

    // It's cheaper to check every populated grid cell on this grid level.
    for (typename map<GridCellIndex<Dimension>, int>::const_iterator gcItr = gridCellPopulations[gridLevel].begin();
         gcItr != gridCellPopulations[gridLevel].end();
         ++gcItr) {
      const GridCellIndex<Dimension> gc = gcItr->first;
      bool useGC = false;
      for (int i = 0; i != Dimension::nDim; ++i) {
        if (abs(gc(i) - gridCell(i)) == gridCellStep) useGC = true;
      }
      if (useGC) result.push_back(gc);
    }

  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Compute the set of grid cells at the given step from a center grid cell.
//------------------------------------------------------------------------------
template<>
inline
set<GridCellIndex<Dim<1> > >
NestedGridRedistributeNodes<Dim<1> >::
computeGridCellRind(const GridCellIndex<Dim<1> >& gridCell,
                    const int gridCellStep) const {
  typedef GridCellIndex<Dim<1> > GridCellType;
  set<GridCellType> result;
  result.insert(GridCellType(gridCell.xIndex() - gridCellStep));
  result.insert(GridCellType(gridCell.xIndex() + gridCellStep));
  ENSURE((result.size() == 1 && gridCellStep == 0) ||
         (result.size() == 2 && gridCellStep > 0));
  return result;
}

template<>
inline
set<GridCellIndex<Dim<2> > >
NestedGridRedistributeNodes<Dim<2> >::
computeGridCellRind(const GridCellIndex<Dim<2> >& gridCell,
                    const int gridCellStep) const {
  REQUIRE(gridCellStep >= 0);
  typedef GridCellIndex<Dim<2> > GridCellType;
  set<GridCellType> result;

  const int gcx = gridCell.xIndex();
  const int gcy = gridCell.yIndex();
  for (int i = 0; i <= gridCellStep; ++i) {
    result.insert(GridCellType(gcx - i, gcy - gridCellStep));
    result.insert(GridCellType(gcx + i, gcy - gridCellStep));
    result.insert(GridCellType(gcx - i, gcy + gridCellStep));
    result.insert(GridCellType(gcx + i, gcy + gridCellStep));

    result.insert(GridCellType(gcx - gridCellStep, gcy - i));
    result.insert(GridCellType(gcx - gridCellStep, gcy + i));
    result.insert(GridCellType(gcx + gridCellStep, gcy - i));
    result.insert(GridCellType(gcx + gridCellStep, gcy + i));
  }

  ENSURE(result.size() == max(1, 8*gridCellStep));
  return result;
}

template<>
inline
set<GridCellIndex<Dim<3> > >
NestedGridRedistributeNodes<Dim<3> >::
computeGridCellRind(const GridCellIndex<Dim<3> >& gridCell,
                    const int gridCellStep) const {
  REQUIRE(gridCellStep >= 0);
  typedef GridCellIndex<Dim<3> > GridCellType;
  set<GridCellType> result;

  const GridCellType gcMin = gridCell - gridCellStep;
  const GridCellType gcMax = gridCell + gridCellStep;

  // Iterate over the dimension we'll hold fixed.
  for (int iDim = 0; iDim != 3; ++iDim) {

    // Choose the independent dimensions we'll be sweeping over.
    int ixDim, iyDim;
    switch(iDim) {
    case 0:
      ixDim = 1;
      iyDim = 2;
      break;
    case 1:
      ixDim = 0;
      iyDim = 2;
      break;
    case 2:
      ixDim = 0;
      iyDim = 1;
      break;
    }

    // Sweep the negative and positive planes simultaneously.
    for (int i = 0; i <= 2*gridCellStep; ++i) {
      for (int j = 0; j <= 2*gridCellStep; ++j) {
        GridCellType gcLow(gcMin);
        GridCellType gcHigh(gcMax);
        gcLow(ixDim) += i;
        gcLow(iyDim) += j;
        gcHigh(ixDim) -= i;
        gcHigh(iyDim) -= j;
        result.insert(gcLow);
        result.insert(gcHigh);
      }
    }

  }

  ENSURE(result.size() == 
         4*(2*gridCellStep)*(2*gridCellStep) + 
         (2*gridCellStep + 1)*(2*gridCellStep + 1) + 
         (2*gridCellStep - 1)*(2*gridCellStep - 1));
  return result;
}
         
//------------------------------------------------------------------------------
// Set the master & coarse neighbor info for each NodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridRedistributeNodes<Dimension>::
setMasterNodeLists(DataBase<Dimension>& dataBase,
                   const GridCellIndex<Dimension>& gridCell,
                   const int gridLevel,
                   std::vector<std::vector<int>>& masterLists,
                   std::vector<std::vector<int>>& coarseNeighbors) const {
  REQUIRE(masterLists.size() == dataBase.numNodeLists());
  REQUIRE(coarseNeighbors.size() == dataBase.numNodeLists());
  unsigned iNodeList = 0;
  for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr, ++iNodeList) {
    NestedGridNeighbor<Dimension>& neighbor = getNestedGridNeighbor(*nodeListItr);
    neighbor.setNestedMasterList(gridCell, gridLevel, masterLists[iNodeList], coarseNeighbors[iNodeList]);
  }
}

//------------------------------------------------------------------------------
// Gather up the unassigned coarse neighbor nodes, filling in the global node 
// indices and work.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NestedGridRedistributeNodes<Dimension>::
gatherAvailableCoarseNodes(const DataBase<Dimension>& dataBase,
                           const vector<DomainNode<Dimension> >& nodeDistribution,
                           const FieldList<Dimension, Scalar>& work,
                           const vector<vector<int>>& localCoarseNeighbors,
                           vector<int>& globalNodeIndices,
                           vector<Scalar>& globalNodeWork) const {

  // Pre-conditions.
  REQUIRE(globalNodeIndices.size() == 0);
  REQUIRE(globalNodeWork.size() == 0);

  // Iterate over the coarse neighbor set.
  for (CoarseNodeIterator<Dimension> nodeItr = dataBase.coarseNodeBegin(localCoarseNeighbors);
       nodeItr != dataBase.coarseNodeEnd();
       ++nodeItr) {

    // Find this node in the node distribution.
    typename vector<DomainNode<Dimension> >::const_iterator domainNodeItr = nodeDistribution.begin();
    while (domainNodeItr != nodeDistribution.end() &&
           (domainNodeItr->nodeListID != nodeItr.fieldID() || domainNodeItr->localNodeID != nodeItr.nodeID())) ++domainNodeItr;
    CHECK(domainNodeItr != nodeDistribution.end());
    CHECK(domainNodeItr->nodeListID == nodeItr.fieldID());
    CHECK(domainNodeItr->localNodeID == nodeItr.nodeID());

    // If this node is still unassigned, add it to the result.
    if (domainNodeItr->domainID == -1) {
      globalNodeIndices.push_back(domainNodeItr->globalNodeID);
      globalNodeWork.push_back(work(nodeItr));
    }
  }

  // Currently every process has built it's own local set of available coarse nodes.
  // Reduce that set to process 0.
  const int procID = this->domainID();
  const int numProcs = this->numDomains();
  if (procID == 0) {

    // Get the other processes info in turn.
    for (int sendProc = 1; sendProc != numProcs; ++sendProc) {
      MPI_Status status1, status2, status3;
      int numRecvNodes;
      MPI_Recv(&numRecvNodes, 1, MPI_INT, sendProc, 120, Communicator::communicator(), &status1);
      vector<int> recvNodeIDs(numRecvNodes);
      MPI_Recv(&(*recvNodeIDs.begin()), numRecvNodes, MPI_INT, sendProc, 121, Communicator::communicator(), &status2);
      vector<double> recvNodeWork(numRecvNodes);
      MPI_Recv(&(*recvNodeWork.begin()), numRecvNodes, MPI_DOUBLE, sendProc, 122, Communicator::communicator(), &status3);

      // Add the other processors nodes to our list.
      globalNodeIndices.reserve(globalNodeIndices.size() + numRecvNodes);
      globalNodeWork.reserve(globalNodeWork.size() + numRecvNodes);
      for (int i = 0; i != numRecvNodes; ++i) {
        CHECK(i < recvNodeIDs.size() && i < recvNodeWork.size());
        globalNodeIndices.push_back(recvNodeIDs[i]);
        globalNodeWork.push_back(recvNodeWork[i]);
      }
    }

  } else {

    // Send our info to the root process.
    CHECK(globalNodeIndices.size() == globalNodeWork.size());
    int numSendNodes = globalNodeIndices.size();
    MPI_Send(&numSendNodes, 1, MPI_INT, 0, 120, Communicator::communicator());
    MPI_Send(&(*globalNodeIndices.begin()), numSendNodes, MPI_INT, 0, 121, Communicator::communicator());
    MPI_Send(&(*globalNodeWork.begin()), numSendNodes, MPI_DOUBLE, 0, 122, Communicator::communicator());

  }
  MPI_Barrier(Communicator::communicator());

  // OK, now the root process has all info, while all other processes have their own
  // local set of coarse nodes only.
  BEGIN_CONTRACT_SCOPE
  ENSURE(globalNodeIndices.size() == globalNodeWork.size());
  sort(globalNodeIndices.begin(), globalNodeIndices.end());
  ENSURE(unique(globalNodeIndices.begin(), globalNodeIndices.end()) == globalNodeIndices.end());
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Assign nodes from the available coarse set to the given domain, until it's
// desired work quota is filled or we exhaust the given set of nodes.
// Returns true if the domain was filled on this pass.  Also update the grid cell
// and grid level populations as nodes are assigned and depleted from these 
// pools.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
NestedGridRedistributeNodes<Dimension>::
assignNodesToDomain(const DataBase<Dimension>& dataBase,
                    const vector<int>& globalNodeIndices,
                    const vector<Scalar>& globalNodeWork,
                    const int currentDomainID,
                    const double targetWork,
                    double& currentWorkSum,
                    vector<DomainNode<Dimension> >& nodeDistribution,
                    vector< map<GridCellIndex<Dimension>, int> >&gridCellPopulations,
                    vector<int>& gridLevelPopulations) const {

  REQUIRE(globalNodeIndices.size() == globalNodeWork.size());
  REQUIRE(currentDomainID >= 0 && currentDomainID < this->numDomains());
  REQUIRE(distinctlyGreaterThan(targetWork, 0.0));
  REQUIRE(currentWorkSum >= 0.0);
  REQUIRE(gridCellPopulations.size() == gridLevelPopulations.size());

  // The usual mpi rank info.
  const int procID = this->domainID();
  const int numProcs = this->numDomains();

  // Iterate until we either exhaust the current pool of coarse nodes
  // or fill the domains work quota.
  const int numGridLevels = gridLevelPopulations.size();
  int numCoarseNodes = globalNodeIndices.size();
  MPI_Bcast(&numCoarseNodes, 1, MPI_INT, 0, Communicator::communicator());
  int i = 0;
  bool domainFull = false;
  while (i < numCoarseNodes && !domainFull) {

    // Only process 0 has the information to decide which node is next.
    int globalNodeID;
    double globalWork;
    if (procID == 0) {
      CHECK(i >= 0 && i < numCoarseNodes);
      globalNodeID = globalNodeIndices[i];
      globalWork = globalNodeWork[i];
//       cerr << "    -> Assigning node " << globalNodeID << " to domain " << currentDomainID << endl;
    }
    MPI_Bcast(&globalNodeID, 1, MPI_INT, 0, Communicator::communicator());
    MPI_Bcast(&globalWork, 1, MPI_DOUBLE, 0, Communicator::communicator());

    // Are we the process that has this node?
    GridCellIndex<Dimension> gc;
    int gridLevel;
    int ownerProc = -1;
    typename vector<DomainNode<Dimension> >::iterator domainNodeItr = nodeDistribution.begin();
    while (domainNodeItr != nodeDistribution.end() &&
           domainNodeItr->globalNodeID != globalNodeID) ++domainNodeItr;
    if (domainNodeItr != nodeDistribution.end()) {

      // Yep, we're the domain that has the node, so assign it.
      ownerProc = procID;
      CHECK(domainNodeItr->domainID == -1);
      domainNodeItr->domainID = currentDomainID;
      
      // Figure out which grid cell/grid level this node is on.
      CHECK(domainNodeItr->nodeListID < dataBase.numNodeLists());
      const NestedGridNeighbor<Dimension>& neighbor = getNestedGridNeighbor(*(dataBase.nodeListBegin() + 
                                                                              domainNodeItr->nodeListID));
      gridLevel = neighbor.gridLevel(domainNodeItr->localNodeID);
      gc = neighbor.nodeInCell()[gridLevel][domainNodeItr->localNodeID];

    }

    // Broadcast the node's grid cell/grid level info to everyone.
    int globalOwnerProc = -1;
    MPI_Allreduce(&ownerProc, &globalOwnerProc, 1, MPI_INT, MPI_MAX, Communicator::communicator());
    CHECK(globalOwnerProc >= 0 && globalOwnerProc < numProcs);
    MPI_Bcast(&gridLevel, 1, MPI_INT, globalOwnerProc, Communicator::communicator());
    CHECK(gridLevel >= 0 && gridLevel < numGridLevels);
    MPI_Bcast(&(gc(0)), Dimension::nDim, MPI_INT, globalOwnerProc, Communicator::communicator());

    // Remove this node from the grid cell/grid level pools of available nodes.
    CHECK(gridLevel < gridCellPopulations.size());
    CHECK(gridCellPopulations[gridLevel].find(gc) != gridCellPopulations[gridLevel].end());
    gridCellPopulations[gridLevel][gc] -= 1;
    gridLevelPopulations[gridLevel] -= 1;

    // Increment the work assigned to this domain, and check if it's full yet.
    currentWorkSum += globalWork;
    if (currentWorkSum >= targetWork) domainFull = true;
    ++i;

  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    // Check that everyone agrees if the domain is full or not.
    int fullFlag = domainFull ? 1 : 0;
    int fullFlagSum;
    MPI_Allreduce(&fullFlag, &fullFlagSum, 1, MPI_INT, MPI_SUM, Communicator::communicator());
    ENSURE(fullFlagSum == numProcs * fullFlag);

    // Check that we're all consistent about the grid level & grid cell populations.
    for (int gl = 0; gl != numGridLevels; ++gl) {
      int glPop = gridLevelPopulations[gl];
      int gcPop = 0;
      for (typename map<GridCellIndex<Dimension>, int>::const_iterator itr = gridCellPopulations[gl].begin();
           itr != gridCellPopulations[gl].end();
           ++itr) gcPop += itr->second;
      CHECK(glPop == gcPop);
      int glPopSum, gcPopSum;
      MPI_Allreduce(&glPop, &glPopSum, 1, MPI_INT, MPI_SUM, Communicator::communicator());
      MPI_Allreduce(&gcPop, &gcPopSum, 1, MPI_INT, MPI_SUM, Communicator::communicator());
      CHECK(glPopSum == numProcs * glPop);
      CHECK(gcPopSum == numProcs * gcPop);
    }
  }
  END_CONTRACT_SCOPE

  // That's it.  Return whether or not the domain has been filled.
  return domainFull;
}

}
