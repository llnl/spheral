//---------------------------------Spheral++----------------------------------//
// VoronoiRedistributeNodes
//
// This algorithm uses the Voronoi tessellation to decide how to domain 
// decompose our points.  The idea is to relax a set of generator points into
// the SPH node distribution -- the generators are attracted to the SPH points
// repelled by one and other.  These generator points then become the seeds to
// draw the Voronoi tessellation about, each cell of which then represents a 
// computational domain.
//
// Created by JMO, Fri Jan 15 09:56:56 PST 2010
//----------------------------------------------------------------------------//
#include "VoronoiRedistributeNodes.hh"
#include "Utilities/DomainNode.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "FieldOperations/binFieldList2Lattice.hh"
#include "NodeList/NodeList.hh"
#include "Kernel/TableKernel.hh"
#include "Kernel/BSplineKernel.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/RedistributionRegistrar.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/boundPointWithinBox.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/PairComparisons.hh"
#include "allReduce.hh"
#include "Communicator.hh"

#include "Utilities/DBC.hh"

#include <algorithm>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <bitset>
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
// Helper method to find the nearest position in a vector of positions to the 
// given one.
//------------------------------------------------------------------------------
template<typename Dimension>
size_t
findNearestGenerator(const typename Dimension::Vector& xi,
                     const vector<typename Dimension::Vector>& generators,
                     const vector<int>& generatorFlags) {

  // This N^2 thing shouldn't be too bad until we get to lots of processors.
  // Then we'll have to do something smarter.
  size_t result = 0;
  double minR2 = DBL_MAX;
  for (size_t igen = 0; igen != generators.size(); ++igen) {
    if (generatorFlags[igen] == 1) {
      const double dr2 = (xi - generators[igen]).magnitude2();
      // integrateThroughMeshAlongSegment<Dimension, double>(workBins, xmin, xmax, ncells, generators[igen], xi);
      if (dr2 < minR2) {
        minR2 = dr2;
        result = igen;
      }
    }
  }
  ENSURE(result < generators.size() or (generators.size() == 0 and result == 0));
  return result;
}

//------------------------------------------------------------------------------
// Compute the center cell position.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Vector
computeCellPosition(const typename Dimension::Vector& xmin,
		    const typename Dimension::Vector& xmax,
		    const unsigned nxCells,
		    const unsigned index);

template<>
inline
Dim<1>::Vector
computeCellPosition<Dim<1> >(const Dim<1>::Vector& xmin,
			     const Dim<1>::Vector& xmax,
			     const unsigned nxCells,
			     const unsigned index) {
  REQUIRE(index < nxCells);
  return Dim<1>::Vector(xmin.x() + (xmax.x() - xmin.x())/nxCells * (index + 0.5));
}

template<>
inline
Dim<2>::Vector
computeCellPosition<Dim<2> >(const Dim<2>::Vector& xmin,
			     const Dim<2>::Vector& xmax,
			     const unsigned nxCells,
			     const unsigned index) {
  REQUIRE(index < nxCells*nxCells);
  const unsigned ix = index % nxCells;
  const unsigned iy = index / nxCells;
  const Dim<2>::Vector xstep = (xmax - xmin)/nxCells;
  return xmin + Dim<2>::Vector(xstep.x() * (ix + 0.5),
                               xstep.y() * (iy + 0.5));
}

template<>
inline
Dim<3>::Vector
computeCellPosition<Dim<3> >(const Dim<3>::Vector& xmin,
			     const Dim<3>::Vector& xmax,
			     const unsigned nxCells,
			     const unsigned index) {
  REQUIRE(index < nxCells*nxCells*nxCells);
  const unsigned iz = index / (nxCells*nxCells);
  const unsigned iy = (index - iz*nxCells*nxCells) / nxCells;
  const unsigned ix = index % nxCells;
  const Dim<3>::Vector xstep = (xmax - xmin)/nxCells;
  return xmin + Dim<3>::Vector(xstep.x() * (ix + 0.5),
                               xstep.y() * (iy + 0.5),
                               xstep.z() * (iz + 0.5));
}

//------------------------------------------------------------------------------
// Compute the cell boundaries.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
computeCellBoundaries(const typename Dimension::Vector& xmin,
		      const typename Dimension::Vector& xmax,
		      const size_t nxCells,
		      const size_t index,
		      typename Dimension::Vector& xcellMin,
		      typename Dimension::Vector& xcellMax);

template<>
inline
void
computeCellBoundaries<Dim<1> >(const Dim<1>::Vector& xmin,
			       const Dim<1>::Vector& xmax,
			       const size_t nxCells,
			       const size_t index,
			       Dim<1>::Vector& xcellMin,
			       Dim<1>::Vector& xcellMax) {
  REQUIRE(index < nxCells);
  typedef Dim<1>::Vector Vector;
  const Vector xstep = (xmax - xmin)/nxCells;
  xcellMin = Vector(xstep.x()*index);
  xcellMax = Vector(xstep.x()*(index + 1));
}

template<>
inline
void
computeCellBoundaries<Dim<2> >(const Dim<2>::Vector& xmin,
			       const Dim<2>::Vector& xmax,
			       const size_t nxCells,
			       const size_t index,
			       Dim<2>::Vector& xcellMin,
			       Dim<2>::Vector& /*xcellMax*/) {
  REQUIRE(index < nxCells*nxCells);
  typedef Dim<2>::Vector Vector;
  const size_t ix = index % nxCells;
  const size_t iy = index / nxCells;
  const Vector xstep = (xmax - xmin)/nxCells;
  xcellMin = Vector(xstep.x()*ix,
		    xstep.y()*iy);
  xcellMin = Vector(xstep.x()*(ix + 1),
		    xstep.y()*(iy + 1));
}

template<>
inline
void
computeCellBoundaries<Dim<3> >(const Dim<3>::Vector& xmin,
			       const Dim<3>::Vector& xmax,
			       const size_t nxCells,
			       const size_t index,
			       Dim<3>::Vector& xcellMin,
			       Dim<3>::Vector& xcellMax) {
  REQUIRE(index < nxCells*nxCells*nxCells);
  typedef Dim<3>::Vector Vector;
  const size_t iz = 2*(index / (nxCells*nxCells));
  const size_t iy = 2*((index - iz/2) / nxCells);
  const size_t ix = 2*(index % nxCells);
  const Vector xstep = (xmax - xmin)/nxCells;
  xcellMin = Vector(xstep.x() * ix,
		    xstep.y() * iy,
		    xstep.z() * iz);
  xcellMax = Vector(xstep.x() * (ix + 1),
		    xstep.y() * (iy + 1),
		    xstep.z() * (iz + 1));
}

//------------------------------------------------------------------------------
// Compute the subcell positions.
//------------------------------------------------------------------------------
// 1-D
inline
vector<Dim<1>::Vector>
computeDaughterPositions(const Dim<1>::Vector& xmin,
                         const Dim<1>::Vector& xmax) {
  typedef Dim<1>::Vector Vector;
  vector<Vector> result;
  const Vector delta = xmax - xmin;
  result.push_back(xmin + 0.25*delta);
  result.push_back(xmin + 0.75*delta);
  ENSURE(result.size() == 2);
  return result;
}

// 2-D
inline
vector<Dim<2>::Vector>
computeDaughterPositions(const Dim<2>::Vector& xmin,
                         const Dim<2>::Vector& xmax) {
  typedef Dim<2>::Vector Vector;
  vector<Vector> result;
  const Vector delta = xmax - xmin;
  result.push_back(Vector(xmin.x() + 0.25*delta.x(), xmin.y() + 0.25*delta.y()));
  result.push_back(Vector(xmin.x() + 0.75*delta.x(), xmin.y() + 0.25*delta.y()));
  result.push_back(Vector(xmin.x() + 0.25*delta.x(), xmin.y() + 0.75*delta.y()));
  result.push_back(Vector(xmin.x() + 0.75*delta.x(), xmin.y() + 0.75*delta.y()));
  ENSURE(result.size() == 4);
  return result;
}

// 3-D
inline
vector<Dim<3>::Vector>
computeDaughterPositions(const Dim<3>::Vector& xmin,
                         const Dim<3>::Vector& xmax) {
  typedef Dim<3>::Vector Vector;
  vector<Vector> result;
  const Vector delta = xmax - xmin;
  result.push_back(Vector(xmin.x() + 0.25*delta.x(), xmin.y() + 0.25*delta.y(), xmin.z() + 0.25*delta.z()));
  result.push_back(Vector(xmin.x() + 0.75*delta.x(), xmin.y() + 0.25*delta.y(), xmin.z() + 0.25*delta.z()));
  result.push_back(Vector(xmin.x() + 0.25*delta.x(), xmin.y() + 0.75*delta.y(), xmin.z() + 0.25*delta.z()));
  result.push_back(Vector(xmin.x() + 0.75*delta.x(), xmin.y() + 0.75*delta.y(), xmin.z() + 0.25*delta.z()));
  result.push_back(Vector(xmin.x() + 0.25*delta.x(), xmin.y() + 0.25*delta.y(), xmin.z() + 0.75*delta.z()));
  result.push_back(Vector(xmin.x() + 0.75*delta.x(), xmin.y() + 0.25*delta.y(), xmin.z() + 0.75*delta.z()));
  result.push_back(Vector(xmin.x() + 0.25*delta.x(), xmin.y() + 0.75*delta.y(), xmin.z() + 0.75*delta.z()));
  result.push_back(Vector(xmin.x() + 0.75*delta.x(), xmin.y() + 0.75*delta.y(), xmin.z() + 0.75*delta.z()));
  ENSURE(result.size() == 8);
  return result;
}

//------------------------------------------------------------------------------
// Compute the node position closest to the cell center.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Vector
computeClosestNodePosition(const typename Dimension::Vector& targetPosition,
                           const vector<DomainNode<Dimension> >& nodes,
                           const int numProcs, 
                           MPI_Comm communicator = Communicator::communicator()) {
  typedef typename Dimension::Vector Vector;

  // First find the local node closest to the center.
  Vector localResult;
  double minr2 = DBL_MAX;
  for (typename vector<DomainNode<Dimension> >::const_iterator itr = nodes.begin();
       itr != nodes.end();
       ++itr) {
    const double r2 = (itr->position - targetPosition).magnitude2();
    if (r2 < minr2) {
      localResult = itr->position;
      minr2 = r2;
    }
  }
  CHECK(minr2 < DBL_MAX);

  // Find the global minimum.
  Vector result;
  minr2 = DBL_MAX;
  for (int sendProc = 0; sendProc != numProcs; ++sendProc) {
    vector<char> buffer;
    packElement(localResult, buffer);
    MPI_Bcast(&buffer.front(), buffer.size(), MPI_CHAR, sendProc, communicator);
    vector<char>::const_iterator itr = buffer.begin();
    Vector xi;
    unpackElement(xi, itr, buffer.end());
    CHECK(itr == buffer.end());
    const double r2 = (xi - targetPosition).magnitude2();
    if (r2 < minr2) {
      result = xi;
      minr2 = r2;
    }
  }
  CHECK(minr2 < DBL_MAX);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
VoronoiRedistributeNodes<Dimension>::
VoronoiRedistributeNodes(double ,
                         const bool workBalance,
                         const bool balanceGenerators,
                         const double tolerance,
                         const unsigned maxIterations):
  RedistributeNodes<Dimension>(),
  mWorkBalance(workBalance),
  mBalanceGenerators(balanceGenerators),
  mTolerance(tolerance),
  mMaxIterations(maxIterations) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
VoronoiRedistributeNodes<Dimension>::
~VoronoiRedistributeNodes() {
}

//------------------------------------------------------------------------------
// The main method of this class.  Call on the Voronoi library to describe an
// optimal partioning of the nodes, and then apply that partitioning.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiRedistributeNodes<Dimension>::
redistributeNodes(DataBase<Dimension>& dataBase,
                  vector<Boundary<Dimension>*> boundaries) {

  // The usual parallel info.
  const int numProcs = this->numDomains();
  const int procID = this->domainID();

  // Get the global IDs.
  const FieldList<Dimension, size_t> globalIDs = globalNodeIDs(dataBase);

  // Compute the work and number density per node.
  const TableKernel<Dimension> W(BSplineKernel<Dimension>(), 100u);
  FieldList<Dimension, Scalar> workField = dataBase.newGlobalFieldList(1.0, "work");
  if (this->workBalance()) { // or mBalanceGenerators) {

    // Enforce boundary conditions for the work computation.
    for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
         nodeListItr != dataBase.nodeListEnd();
         ++nodeListItr) {
      (*nodeListItr)->numGhostNodes(0);
      (*nodeListItr)->neighbor().updateNodes();
    }
    for (typename vector<Boundary<Dimension>*>::iterator boundaryItr = boundaries.begin(); 
         boundaryItr != boundaries.end();
         ++boundaryItr) {
      (*boundaryItr)->setAllGhostNodes(dataBase);
      (*boundaryItr)->finalizeGhostBoundary();
      for (typename DataBase<Dimension>::FluidNodeListIterator nodeListItr = dataBase.fluidNodeListBegin();
           nodeListItr != dataBase.fluidNodeListEnd(); 
           ++nodeListItr) (*nodeListItr)->neighbor().updateNodes();
    }

    // Update the connectivity.
    dataBase.updateConnectivityMap(false, false, false);

    // Get the local description of the domain distribution, with the work per node filled in.
    if (this->workBalance()) workField = this->workPerNode(dataBase, 1.0);
  }

  // Print the beginning statistics.
  std::string stats0 = this->gatherDomainDistributionStatistics(workField);
  if (procID == 0) cout << "VoronoiRedistributeNodes: INITIAL node distribution statistics:" << endl
                        << stats0 << endl;

  // Now we can get the node distribution description.
  vector<DomainNode<Dimension> > nodeDistribution = this->currentDomainDecomposition(dataBase, globalIDs, workField);
  const size_t numNodes = nodeDistribution.size();
  const size_t numNodesGlobal = allReduce((uint64_t) numNodes, SPHERAL_OP_SUM);
  const size_t avgNumNodes = numNodesGlobal/numProcs;
  CHECK(numNodes > 0);

  // Clear out any ghost nodes.
  // We won't bother to update the neighbor info at this point -- we don't need 
  // it for this algorithm, so we just update it when we're done.
  for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) (*nodeListItr)->numGhostNodes(0);

  // Get the bounding boxes for the node set.
  Vector xmin, xmax;
  dataBase.boundingBox(xmin, xmax);

  // Define the the length scale we use to determine when the generator positions have converged.
  const double tol = (xmax - xmin).minElement() * mTolerance;
  if (procID == 0) cout << "VoronoiRedistributeNodes: Found bounding box of " << xmin << " " << xmax << endl
                        << "                          yielding generator convergence tolerance of " << tol << endl;

  // Determine the average work per generator.
  const Scalar totWork = workField.sumElements();
  const Scalar avgWork = totWork/numProcs;
  
  // Try and distribute the seeds according to where the work is by an AMR like zooming
  // in of the work distribution in bins.
  vector<Vector> generators(numProcs, 0.5*(xmin + xmax));
  vector<pair<Vector, Vector> > generatorBounds(numProcs, make_pair(xmin, xmax));
  const size_t numDaughters = (1U << Dimension::nDim);
  vector<unsigned> ncells(Dimension::nDim, 2U);
  {
    // Everyone starts out in the same bin.
    vector<vector<size_t> > generatorsInParents(1);
    vector<pair<Vector, Vector> > parentCells;
    for (size_t igen = 0; (int)igen != numProcs; ++igen) generatorsInParents[0].push_back(igen);
    parentCells.push_back(make_pair(xmin, xmax));

    // Descend until we either get each generator in an individual cell or hit the maximum number
    // of allowed levels
    size_t level = 0;
    size_t numRemainingGenerators = numProcs;
    while (generatorsInParents.size() > 0 and level < mMaxIterations) {
      CHECK(generatorsInParents.size() == parentCells.size());
      ++level;
      numRemainingGenerators = 0;
      vector<vector<size_t> > newGeneratorsInParents;
      vector<pair<Vector, Vector> > newParentCells;

      // Walk each parent cell.
      for (size_t k = 0; k != generatorsInParents.size(); ++k) {
        const vector<size_t>& gens = generatorsInParents[k];
        const size_t numParentGenerators = gens.size();
        const Vector& xminParent = parentCells[k].first;
        const Vector& xmaxParent = parentCells[k].second;

        // Break up the work in this parent cell into 2^ndim subcells.
        const vector<double> workBins = binFieldList2Lattice(workField, xminParent, xmaxParent, ncells);
        const double workParent = accumulate(workBins.begin(), workBins.end(), 0.0);
        CHECK(workBins.size() == numDaughters);
        CHECK(workParent > 0.0);
        
        // Sort the daughters by work.
        vector<pair<double, size_t> > sortedDaughters;
        for (size_t k = 0; k != numDaughters; ++k) sortedDaughters.push_back(make_pair(workBins[k], k));
        sort(sortedDaughters.begin(), sortedDaughters.end(), ComparePairsByFirstElementInDecreasingOrder<pair<double, size_t> >());

        // Figure out how many generators we're assigning to each daughter.
        vector<size_t> numGensForDaughter(numDaughters, 0);
        {
          size_t numAssigned = 0;
          size_t k = 0;
          while (numAssigned < numParentGenerators and k < numDaughters and sortedDaughters[k].first > 0.0) {
            const double work = sortedDaughters[k].first;
            const size_t kdaughter = sortedDaughters[k].second;
            numGensForDaughter[kdaughter] = min(numParentGenerators - numAssigned, size_t(work*safeInv(workParent, 1.0e-30)*numParentGenerators + 1));
            numAssigned += numGensForDaughter[kdaughter];
            ++k;
          }

          // Pick up any stragglers.
          VERIFY(numAssigned == numParentGenerators);
          if (numAssigned < numParentGenerators) {
            k = 0;
            while (numAssigned < numParentGenerators and k < numDaughters and sortedDaughters[k].first > 0.0) {
              const size_t kdaughter = sortedDaughters[k].second;
              ++numGensForDaughter[kdaughter];
              ++numAssigned;
              ++k;
            }
          }
          VERIFY(numAssigned == numParentGenerators);
        }

//         // Iterate until we have assigned each of the parent generators to a daughter.
//         vector<size_t> numGensForDaughter(numDaughters, 0);
//         int numRemaining = numParentGenerators;
//         int iteration = 0;
//         while (iteration < mMaxIterations and numRemaining > 0) {
//           ++iteration;
//           int numStillRemaining = numRemaining;
//           for (size_t kdaughter = 0; kdaughter != workBins.size(); ++kdaughter) {
//             const size_t delta = min(size_t(numStillRemaining), size_t(workBins[kdaughter]*numRemaining*safeInv(workParent, 1.0e-30) + 0.5 + 0.5*double(iteration)/double(max(1U, mMaxIterations - 1))));
//             numGensForDaughter[kdaughter] += delta;
//             numStillRemaining -= delta;
//             CHECK(numStillRemaining >= 0);
//           }
//           numRemaining = numStillRemaining;
//         }
//         VERIFY(numRemaining == 0);

        // Find the positions of the daughter cells.
        const vector<Vector> daughterPositions = computeDaughterPositions(xminParent, xmaxParent);
        CHECK(daughterPositions.size() == numDaughters);

        // Now assign the generator positions and the next generation of parents.
        vector<size_t>::const_iterator genItr = gens.begin();
        for (size_t kdaughter = 0; kdaughter != numDaughters; ++kdaughter) {
          CHECK(genItr + numGensForDaughter[kdaughter] <= gens.end());
          const Vector dcell = 0.25*(xmaxParent - xminParent);
          const Vector xminDaughter = daughterPositions[kdaughter] - dcell;
          const Vector xmaxDaughter = daughterPositions[kdaughter] + dcell;
          if (numGensForDaughter[kdaughter] == 1) {
            generators[*genItr] = computeClosestNodePosition<Dimension>(0.5*(xminDaughter + xmaxDaughter),
                                                                        nodeDistribution, numProcs);
            generatorBounds[*genItr] = make_pair(xminDaughter, xmaxDaughter);
            CHECK(testPointInBox(generators[*genItr], xminDaughter, xmaxDaughter));
          } else if (numGensForDaughter[kdaughter] > 1) {
            for (vector<size_t>::const_iterator itr = genItr; itr != genItr + numGensForDaughter[kdaughter]; ++itr) {
              generators[*itr] = daughterPositions[kdaughter];
              generatorBounds[*itr] = make_pair(xminDaughter, xmaxDaughter);
            }
            newGeneratorsInParents.push_back(vector<size_t>(genItr, genItr + numGensForDaughter[kdaughter]));
            newParentCells.push_back(make_pair(xminDaughter, xmaxDaughter));
            numRemainingGenerators += numGensForDaughter[kdaughter];
          }
          genItr += numGensForDaughter[kdaughter];
        }
        CHECK(genItr == gens.end());
      }

      // Assign the next generation.
      CHECK(newGeneratorsInParents.size() == newParentCells.size());
      generatorsInParents = newGeneratorsInParents;
      parentCells = newParentCells;
      if (procID == 0) cout << "   Generation " << level << " : "
                            << numRemainingGenerators << " generators remaining in " 
                            << generatorsInParents.size() << " cells."
                            << endl;
    }
    VERIFY(numRemainingGenerators == 0);

//     // Are there still remaining degeneracies in the generator positions?
//     if (numRemainingGenerators > 0) {
//       if (procID == 0) cout << "  --> Breaking up " << numRemainingGenerators 
//                             << " degeneracies in intial generator positions."
//                             << endl;
//       for (vector<vector<size_t> >::const_iterator cellItr = generatorsInParents.begin();
//            cellItr != generatorsInParents.end();
//            ++cellItr) {
//         for (vector<size_t>::const_iterator genItr = cellItr->begin();
//              genItr != cellItr->end();
//              ++genItr) {
//           const size_t igen = *genItr;
//           if (procID == igen) generators[igen] = nodeDistribution[numNodes/2].position;
//           vector<char> buffer;
//           packElement(generators[igen], buffer);
//           MPI_Bcast(&buffer.front(), buffer.size(), MPI_CHAR, igen, Communicator::communicator());
//           vector<char>::const_iterator itr = buffer.begin();
//           unpackElement(generators[igen], itr, buffer.end());
//           CHECK(itr == buffer.end());
//         }
//       }
//     }

  }

  // Copy the initial generator distribution.
  const vector<Vector> startingGenerators(generators);

//   // Stage 1:  Lloyds algorithm iteration of the generator positions.
//   // Choose the initial positions of the generators randomly, one per domain.
//   vector<Vector> generators(numProcs);
//   for (size_t igen = 0; igen != numProcs; ++igen) {
//     if (procID == igen) generators[igen] = nodeDistribution[numNodes/2].position;
//     vector<char> buffer;
//     packElement(generators[igen], buffer);
//     MPI_Bcast(&buffer.front(), buffer.size(), MPI_CHAR, igen, Communicator::communicator());
//     Vector xgen;
//     vector<char>::const_iterator itr = buffer.begin();
//     unpackElement(xgen, itr, buffer.end());
//     CHECK(itr == buffer.end());
//     generators[igen] = xgen;
//   }

  // Set all nodes as unassigned.
  for (size_t i = 0; i != nodeDistribution.size(); ++i) nodeDistribution[i].domainID = -1;

  // Assign the Voronoi distribution based on the initial seeds.
  vector<double> generatorWork(generators.size(), 0.0);
  vector<int> generatorFlags(generators.size(), 1);
  double minWork, maxWork;
  unsigned minNodes, maxNodes;
  assignNodesToGenerators(generators,
                          generatorFlags,
                          generatorWork,
                          nodeDistribution,
                          minWork,
                          maxWork,
                          minNodes,
                          maxNodes);
                          

  // Now iterate the generators until we either converge or hit the max iterations.
  size_t iteration = 0;
  double maxDeltaGenerator = 10.0*tol;
  double oldWorkRatio = maxWork*safeInv(minWork);
  double workRatio = 10.0*oldWorkRatio;
  while (iteration < mMaxIterations and
         maxDeltaGenerator > tol and
         abs(workRatio*safeInv(oldWorkRatio) - 1.0) > 0.001) {

    // Remember the starting generators.
    const vector<Vector> generators0(generators);

    // Iterate until we either have the desired work distribution for this pass or
    // all generators have been flagged.
    generatorFlags = vector<int>(generators.size(), 1);
    while (minWork/maxWork < 1.1 and accumulate(generatorFlags.begin(), generatorFlags.end(), 0) > 1) {

      // Cull out nodes from generators that have too much work.
      cullGeneratorNodesByWork(generators, generatorWork, avgWork, generatorFlags, nodeDistribution);

      // Reassign the loose nodes.
      assignNodesToGenerators(generators, generatorFlags, generatorWork, nodeDistribution,
                              minWork, maxWork, minNodes, maxNodes);

    }

    // Determine the new generators.
    computeCentroids(nodeDistribution, generators);

    // Reestablish the proper Voronoi distribution based on this iterations
    // generators.
    generatorWork = vector<double>(generators.size(), 0.0);
    generatorFlags = vector<int>(generators.size(), 1);
    for (size_t i = 0; i != nodeDistribution.size(); ++i) nodeDistribution[i].domainID = -1;
    assignNodesToGenerators(generators, generatorFlags, generatorWork, nodeDistribution, 
                            minWork, maxWork, minNodes, maxNodes);

    // How much did the generators shift?
    maxDeltaGenerator = 0.0;
    for (size_t igen = 0; igen != generators.size(); ++igen) {
      maxDeltaGenerator = max(maxDeltaGenerator, (generators[igen] - generators0[igen]).magnitude2());
    }
    maxDeltaGenerator = sqrt(maxDeltaGenerator);

    // How much did the work distribution change?
    oldWorkRatio = workRatio;
    workRatio = maxWork*safeInv(minWork);

    // Report this iterations statistics.
    if (procID == 0) cout << "VoronoiRedistributeNodes: Lloyds iteration " << iteration << endl
                          << "                          max change:  " << maxDeltaGenerator << endl
                          << "                          work ratio change:  " << workRatio << " " << oldWorkRatio << " " << abs(workRatio*safeInv(oldWorkRatio) - 1.0) << endl
                          << "                          [min, max, avg] work      [" << minWork << ", " << maxWork << ", " << avgWork << "]" << endl
                          << "                          [min, max, avg] num nodes [" << minNodes << ", " << maxNodes << ", " << avgNumNodes << "]" << endl;
    if (minWork == 0.0) {
      if (procID == 0) {
        cerr << "ERROR:  zero work associated with the following generators:" << endl;
        for (size_t k = 0; (int)k != numProcs; ++k) {
          if (generatorWork[k] == 0.0) {
            cerr << "    ----->  " << generators[k] << endl;
            generators[k] = startingGenerators[k];
          }
        }
      }
    }
    //VERIFY(minWork > 0);
    ++iteration;
  }

  // Redistribute nodes between domains.
  CHECK(this->validDomainDecomposition(nodeDistribution, dataBase));
  this->enforceDomainDecomposition(nodeDistribution, dataBase);

  // Reinitialize neighbor info.
  for (typename DataBase<Dimension>::NodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    (*nodeListItr)->neighbor().updateNodes();
  }

  // Notify everyone that the nodes have just been shuffled around.
  RedistributionRegistrar::instance().broadcastRedistributionNotifications();

  // Print the final statistics.
  std::string stats1 = this->gatherDomainDistributionStatistics(workField);
  if (Process::getRank() == 0) cout << "VoronoiRedistributeNodes: FINAL node distribution statistics:" << endl
                                    << stats1 << endl;
}

//------------------------------------------------------------------------------
// Split the given node distribution into ncells sets using a Voronoi Lloyds type
// algorithm.  The result is returned by setting the domainID attribute of each
// DomainNode to a value in the range [0, ncells[.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiRedistributeNodes<Dimension>::
computeCentroids(const vector<DomainNode<Dimension> >& nodes,
                 vector<typename Dimension::Vector>& generators) const {

  const int numProcs = this->numDomains();
  const size_t numGenerators = generators.size();
  REQUIRE((int)numGenerators == numProcs);

  // Initializations.
  const vector<Vector> generators0(generators);
  generators = vector<Vector>(numGenerators, Vector::zero);

  // Iterate over the nodes, assigning each to it's nearest generator.
  vector<double> normalization(numGenerators, 0.0);
  for (typename vector<DomainNode<Dimension> >::const_iterator itr = nodes.begin();
       itr != nodes.end();
       ++itr) {
    const DomainNode<Dimension>& node = *itr;
    const size_t igen = node.domainID;
    CHECK(igen < numGenerators);
    const double thpt = node.work;
    generators[igen] += thpt*node.position;
    normalization[igen] += thpt;
  }

  // Reduce the generator info across processors.
  vector<char> localBuffer;
  for (size_t igen = 0; igen != numGenerators; ++igen) {
    packElement(generators[igen], localBuffer);
    packElement(normalization[igen], localBuffer);
    generators[igen] = Vector::zero;
    normalization[igen] = 0.0;
  }
  for (size_t sendProc = 0; (int)sendProc != numProcs; ++sendProc) {
    vector<char> buffer = localBuffer;
    MPI_Bcast(&buffer.front(), buffer.size(), MPI_CHAR, sendProc, Communicator::communicator());
    vector<char>::const_iterator itr = buffer.begin();
    for (size_t igen = 0; igen != numGenerators; ++igen) {
      Vector ri;
      double normi;
      unpackElement(ri, itr, buffer.end());
      unpackElement(normi, itr, buffer.end());
      generators[igen] += ri;
      normalization[igen] += normi;
    }
    CHECK(itr == buffer.end());
  }

  // Normalize the new generator positions and compute our final statistics.
  // We also force the generator to the position of the nearest node, ensuring
  // that all generators are actually in the the node distribution and have 
  // at least some work.
  for (size_t igen = 0; igen != numGenerators; ++igen) {
    generators[igen] = 0.25*generators[igen]*safeInv(normalization[igen]) + 0.75*generators0[igen];
    generators[igen] = computeClosestNodePosition<Dimension>(generators[igen], nodes, numProcs);
  }
}

//------------------------------------------------------------------------------
// Split the given node distribution into ncells sets using a Voronoi algorithm.
// The result is returned by setting the domainID attribute of each
// DomainNode to a value in the range [0, ncells).
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiRedistributeNodes<Dimension>::
assignNodesToGenerators(const vector<typename Dimension::Vector>& generators,
                        const vector<int>& generatorFlags,
                        vector<double>& generatorWork,
                        vector<DomainNode<Dimension> >& nodes,
                        double& minWork,
                        double& maxWork,
                        unsigned& minNodes,
                        unsigned& maxNodes) const {

  const int numProcs = this->numDomains();

  // Initializations.
  minWork = DBL_MAX;
  maxWork = -DBL_MAX;
  minNodes = INT_MAX;
  maxNodes = 0;
  const size_t numGenerators = generators.size();
  generatorWork = vector<double>(numGenerators, 0.0);
  vector<unsigned> numNodesPerGenerator(numGenerators, 0);

  // Iterate over the nodes, assigning unassigned nodes to the nearest available generator.
  for (typename vector<DomainNode<Dimension> >::iterator itr = nodes.begin();
       itr != nodes.end();
       ++itr) {
    DomainNode<Dimension>& node = *itr;
    if (node.domainID == -1) {
      const size_t igen = findNearestGenerator<Dimension>(node.position, generators, generatorFlags);
      node.domainID = igen;
    }
    const size_t igen = node.domainID;
    CHECK(igen < numGenerators);
    generatorWork[igen] += node.work;
    ++numNodesPerGenerator[igen];
  }

  // Reduce the generator info across processors.
  vector<char> localBuffer;
  for (size_t igen = 0; igen != numGenerators; ++igen) {
    packElement(generatorWork[igen], localBuffer);
    packElement(numNodesPerGenerator[igen], localBuffer);
    generatorWork[igen] = 0.0;
    numNodesPerGenerator[igen] = 0;
  }
  for (size_t sendProc = 0; (int)sendProc != numProcs; ++sendProc) {
    vector<char> buffer = localBuffer;
    MPI_Bcast(&buffer.front(), buffer.size(), MPI_CHAR, sendProc, Communicator::communicator());
    vector<char>::const_iterator itr = buffer.begin();
    for (size_t igen = 0; igen != numGenerators; ++igen) {
      double worki;
      unsigned ni;
      unpackElement(worki, itr, buffer.end());
      unpackElement(ni, itr, buffer.end());
      generatorWork[igen] += worki;
      numNodesPerGenerator[igen] += ni;
    }
    CHECK(itr == buffer.end());
  }

  // Compute our final statistics.
  for (size_t igen = 0; igen != numGenerators; ++igen) {
    minWork = min(minWork, generatorWork[igen]);
    maxWork = max(maxWork, generatorWork[igen]);
    minNodes = min(minNodes, numNodesPerGenerator[igen]);
    maxNodes = max(maxNodes, numNodesPerGenerator[igen]);
  }
}

//------------------------------------------------------------------------------
// Cull or unassign nodes from generators with too much work.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiRedistributeNodes<Dimension>::
cullGeneratorNodesByWork(const vector<typename Dimension::Vector>& generators,
                         const vector<double>& generatorWork,
                         const double targetWork,
                         vector<int>& generatorFlags,
                         vector<DomainNode<Dimension> >& nodes) const {

  // Pre-conditions.
  const size_t numGenerators = generators.size();
  REQUIRE(targetWork > 0.0);
  REQUIRE(generatorWork.size() == numGenerators);

  // How many generators are we starting with?
  const int numOpenGenerators = accumulate(generatorFlags.begin(), generatorFlags.end(), 0);

  // Walk the generators.
  for (size_t igen = 0; igen != numGenerators; ++igen) {

    // Do we need to examine this generator?
    if (generatorFlags[igen] == 1 and 
        generatorWork[igen] >= targetWork) {
      generatorFlags[igen] = 0;

      // Try to cull out nodes from this generator so we don't overshoot the target work.
      // Find the nodes associated with this generator on this processor.
      // Sort them in increasing distance from the generator.
      typedef pair<size_t, double> PairType;
      vector<PairType> distances;
      for (size_t i = 0; i != nodes.size(); ++i) {
        if (nodes[i].domainID == (int)igen) distances.push_back(make_pair(i, (nodes[i].position - generators[igen]).magnitude2()));
      }
      sort(distances.begin(), distances.end(), ComparePairsBySecondElement<PairType>());

      // Find the global range of distances from the generator.
      double rmin = allReduce((distances.size() > 0 ? distances.front().second : DBL_MAX), SPHERAL_OP_MIN);
      double rmax = allReduce((distances.size() > 0 ? distances.back().second  : 0.0),     SPHERAL_OP_MAX);

      // Bisect for the appropriate radius to reject nodes.
      const double worktol = max(1.0e-10, 0.01*targetWork);
      const double rtol = 1.0e-10*max(1.0, rmax - rmin);
      double rreject = rmax;
      double currentWork = generatorWork[igen];
      while ((abs(currentWork - targetWork) >  worktol) and ((rmax - rmin) > rtol)) {
        rreject = 0.5*(rmin + rmax);
        double localWork = 0.0;
        vector<PairType>::const_iterator itr = distances.begin();
        while (itr != distances.end() and itr->second < rreject) {
          localWork += nodes[itr->first].work;
          ++itr;
        }
        currentWork = allReduce(localWork, SPHERAL_OP_SUM);
        if (currentWork < targetWork) {
          rmin = rreject;
        } else {
          rmax = rreject;
        }
      }

      // Now go through and unassign any nodes from this generator that are outside the
      // rejection threshold.
      typename vector<PairType>::iterator lowerItr = lower_bound(distances.begin(), 
                                                                 distances.end(),
                                                                 rreject,
                                                                 ComparePairsBySecondElement<PairType>());
      for (typename vector<PairType>::iterator itr = lowerItr;
           itr != distances.end();
           ++itr) nodes[itr->first].domainID = -1;
    }
  }

  // Did we actually assign any generators this pass?  If not, we're stuck!
  // Flag all generators as done.
  if (accumulate(generatorFlags.begin(), generatorFlags.end(), 0) == numOpenGenerators) {
    generatorFlags = vector<int>(numGenerators, 0);
  }
}

//------------------------------------------------------------------------------
// Find the Voronoi cells adjacent to the specified one.
// If the given cell is one of the vertices, reject it.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<size_t>
VoronoiRedistributeNodes<Dimension>::
findNeighborGenerators(const size_t igen,
                       const vector<Vector>& generators) const {

  // Are there enough generators to construct a meaningful volume?
  vector<size_t> result;
  if (generators.size() < (Dimension::nDim + 1)) {
    for (size_t i = 0; i != generators.size(); ++i) {
      if (i != igen) result.push_back(i);
    }
    return result;
  }

  // Construct the convex hull of the inverse distance to the other generators.
  vector<Vector> inverseDistance;
  for (size_t jgen = 0; jgen != generators.size(); ++jgen) {
    if (jgen != igen) {
      const Vector delta = generators[jgen] - generators[igen];
      const Scalar rmag2 = delta.magnitude2();
      CHECK(rmag2 > 0.0);
      inverseDistance.push_back(delta * safeInv(rmag2, 1.0e-4));
    } else {
      inverseDistance.push_back(Vector::zero);
    }
  }
  typedef typename Dimension::ConvexHull ConvexHull;
  ConvexHull invHull(inverseDistance);
  const vector<Vector> hullVertices = invHull.vertices();

  // Select the generators that correspond to the vertices of the inverse hull.
  for (size_t jgen = 0; jgen != generators.size(); ++jgen) {
    if (jgen != igen) {
      for (typename vector<Vector>::const_iterator itr = hullVertices.begin();
           itr != hullVertices.end();
           ++itr) {
        if (fuzzyEqual((inverseDistance[jgen] - *itr).magnitude2(), 0.0, 1.0e-10)) result.push_back(jgen);
      }
    }
  }

  return result;
}

//------------------------------------------------------------------------------
// Flag for whether we should compute the work per node or strictly balance by
// node count.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
VoronoiRedistributeNodes<Dimension>::
workBalance() const {
  return mWorkBalance;
}

template<typename Dimension>
void
VoronoiRedistributeNodes<Dimension>::
workBalance(bool val) {
  mWorkBalance = val;
}

//------------------------------------------------------------------------------
// Flag to determine if we try to work balance the generators as part of the 
// iteration.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
VoronoiRedistributeNodes<Dimension>::
balanceGenerators() const {
  return mBalanceGenerators;
}

template<typename Dimension>
void
VoronoiRedistributeNodes<Dimension>::
balanceGenerators(bool val) {
  mBalanceGenerators = val;
}

//------------------------------------------------------------------------------
// Tolerance used to determine convergence of the generators.
//------------------------------------------------------------------------------
template<typename Dimension>
double
VoronoiRedistributeNodes<Dimension>::
tolerance() const {
  return mTolerance;
}

template<typename Dimension>
void
VoronoiRedistributeNodes<Dimension>::
tolerance(double val) {
  mTolerance = val;
}

//------------------------------------------------------------------------------
// Maximum allowed iterations to try and converge the generator positions.
//------------------------------------------------------------------------------
template<typename Dimension>
unsigned
VoronoiRedistributeNodes<Dimension>::
maxIterations() const {
  return mMaxIterations;
}

template<typename Dimension>
void
VoronoiRedistributeNodes<Dimension>::
maxIterations(unsigned val) {
  mMaxIterations = val;
}

}

