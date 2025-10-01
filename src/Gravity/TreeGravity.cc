//---------------------------------Spheral++----------------------------------//
// TreeGravity -- An implementation of the tree n-body gravity solver.
// Based on the original 3D only OctTreeGravity.
//
// Created by JMO, 2013-06-12
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "TreeGravity.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/packElement.hh"
#include "Distributed/allReduce.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/PairComparisons.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Distributed/Communicator.hh"
#include "Utilities/DBC.hh"
#include "Utilities/Hashes.hh"

#include "config.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include <algorithm>
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

namespace {
//------------------------------------------------------------------------------
// Encapculate some specialized dimension specific things we need for this class.
//------------------------------------------------------------------------------
template<typename Dimension> struct TreeDimensionTraits;

//..............................................................................
// 2D
//..............................................................................
template<>
struct TreeDimensionTraits<Dim<2> > {
  typedef Dim<2>::Vector Vector;

  static Vector cellCenter(const Vector& xmin,
                           const unsigned ix, 
                           const unsigned iy, 
                           const unsigned /*iz*/,
                           const double cellsize) {
    return Vector(xmin.x() + (ix + 0.5)*cellsize, 
                  xmin.y() + (iy + 0.5)*cellsize);
  }

  static double forceLaw(const double r2) { return 1.0/sqrt(r2); }
  static double potentialLaw(const double r2) { return log(sqrt(r2)); }
};

//..............................................................................
// 3D
//..............................................................................
template<>
struct TreeDimensionTraits<Dim<3> > {
  typedef Dim<3>::Vector Vector;

  static Vector cellCenter(const Vector& xmin,
                           const unsigned ix, 
                           const unsigned iy, 
                           const unsigned iz,
                           const double cellsize) {
    return Vector(xmin.x() + (ix + 0.5)*cellsize,
                  xmin.y() + (iy + 0.5)*cellsize,
                  xmin.z() + (iz + 0.5)*cellsize);
  }

  static double forceLaw(const double r2) { return 1.0/r2; }
  static double potentialLaw(const double r2) { return 1.0/sqrt(r2); }
};

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TreeGravity<Dimension>::
TreeGravity(const double G,
            const double softeningLength,
            const double opening,
            const double ftimestep,
            const GravityTimeStepType timeStepChoice):
  GenericBodyForce<Dimension >(),
  mG(G),
  mSofteningLength(softeningLength),
  mOpening2(opening*opening),
  mftimestep(ftimestep),
  mBoxLength(0.0),
  mTimeStepChoice(timeStepChoice),
  mXmin(),
  mXmax(),
  mTree(),
  mPotential(FieldStorageType::CopyFields),
  mExtraEnergy(0.0),
  mNodeListMax(0),
  mimax(0),
  mDtMinAcc(0.0),
  mRhoMax(0.0),
  mRestart(registerWithRestart(*this)) {
  VERIFY(G > 0.0);
  VERIFY(opening > 0.0);
  VERIFY(softeningLength > 0.0);
  VERIFY(ftimestep > 0.0);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TreeGravity<Dimension>::
~TreeGravity() {
}

//------------------------------------------------------------------------------
// Register some extra state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeGravity<Dimension>::
registerState(DataBase<Dimension >& dataBase,
              State<Dimension >& state) {
  GenericBodyForce<Dimension >::registerState(dataBase, state);
  state.enroll(mPotential);
}

//------------------------------------------------------------------------------
// Evaluate the forces and return our time derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
TreeGravity<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension >& state,
                    StateDerivatives<Dimension >& derivs) const {
  const double softLength2 = mSofteningLength*mSofteningLength;

  // Access the pertinent fields in the database.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const size_t numNodeLists = position.numFields();

  // Get the acceleration and position change vectors we'll be modifying.
  FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Vector> DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);

  // Zero out the total gravitational potential energy.
  mPotential = 0.0;

  // Initialize the time-step voting data.
  FieldList<Dimension, std::vector<Scalar> > interactionMasses = dataBase.newGlobalFieldList(vector<Scalar>(), "gravity dt interaction masses");
  FieldList<Dimension, std::vector<Vector> > interactionPositions = dataBase.newGlobalFieldList(vector<Vector>(), "gravity dt interaction positions");
  FieldList<Dimension, std::pair<LevelKey, CellKey> > homeBuckets = dataBase.newGlobalFieldList(pair<LevelKey, CellKey>(), "gravity dt home buckets");
  mDtMinAcc = std::numeric_limits<Scalar>::max();

  // Prepare the flags to remember which cells have terminated for each node.
  CompletedCellSet cellsCompleted;
  for (unsigned nodeListi = 0; nodeListi != mass.numFields(); ++nodeListi) {
    for (unsigned i = 0; i != mass[nodeListi]->numInternalElements(); ++i) {
      cellsCompleted[NodeID(nodeListi, i)] = vector<std::unordered_set<CellKey> >(num1dbits);
    }
  }

#ifdef ENABLE_MPI

  // Get the processor information.
  const unsigned rank = Process::getRank();
  const unsigned numProcs = Process::getTotalNumberOfProcesses();

  // Pack up the local tree.
  vector<char> localBuffer, buffer;
  this->serialize(mTree, localBuffer);
  unsigned localBufSize = localBuffer.size();

  // Launch our sends to all other processors.  This may be a bit aggressive.... :)
  vector<MPI_Request> sendRequests;
  sendRequests.reserve(2*numProcs);
  for (unsigned otherProc = 0; otherProc != numProcs; ++otherProc) {
    if (otherProc != rank) {
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&localBufSize, 1, MPI_UNSIGNED, otherProc, 1, Communicator::communicator(), &sendRequests.back());
      if (localBufSize > 0) {
        sendRequests.push_back(MPI_Request());
        MPI_Isend(&localBuffer.front(), localBufSize, MPI_CHAR, otherProc, 2, Communicator::communicator(), &sendRequests.back());
      }
    }
  }
  CHECK(sendRequests.size() <= 2*numProcs);

  // Now walk the other processes and get their trees to add to our own.
  unsigned bufSize;
  MPI_Status recvStatus;
  vector<char>::const_iterator bufItr;
  Tree tree;
  for (unsigned otherProc = 0; otherProc != numProcs; ++otherProc) {
    if (otherProc != rank) {
      MPI_Recv(&bufSize, 1, MPI_UNSIGNED, otherProc, 1, Communicator::communicator(), &recvStatus);
      if (bufSize > 0) {
        buffer = vector<char>(bufSize);
        MPI_Recv(&buffer.front(), bufSize, MPI_CHAR, otherProc, 2, Communicator::communicator(), &recvStatus);
        tree = Tree();
        bufItr = buffer.begin();
        this->deserialize(tree, bufItr, buffer.end());
        CHECK(bufItr == buffer.end());
        mDtMinAcc = min(mDtMinAcc, 
                        applyTreeForces(tree, mass, position, DvDt, mPotential, interactionMasses, interactionPositions, homeBuckets, cellsCompleted));
      }
    }
  }

#endif

  // Apply the forces from our local tree.
  mDtMinAcc = min(mDtMinAcc,
                  applyTreeForces(mTree, mass, position, DvDt, mPotential, interactionMasses, interactionPositions, homeBuckets, cellsCompleted));

  // If we're using the dynamical time step choice, make the calculations for the next
  // time step now.
  // The following is an implementation (perhaps not exact) of the algorithm 
  // described in Zemp et al., MNRAS, 376, 273-286 (2007).
  // The idea is to come up with the most restrictive dynamical time per particle
  // by find the dominant potentials for each point.
  if (mTimeStepChoice == GravityTimeStepType::DynamicalTime) {
    mRhoMax = 0.0;
    mNodeListMax = -1;
    mimax = -1;
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const size_t n = position[nodeListi]->numInternalElements();
      for (unsigned i = 0; i != n; ++i) {
        const Vector& xi = position(nodeListi, i);
        const vector<Scalar>& mbucketsi = interactionMasses(nodeListi, i);
        const vector<Vector>& xbucketsi = interactionPositions(nodeListi, i);
        const unsigned nbuckets = mbucketsi.size();
        CHECK(nbuckets > 0);
        CHECK(xbucketsi.size() == nbuckets);

        // First compute the full set of effective densities for this point, and sort them
        // in decreasing order to pick out the top 0.5%.
        vector<pair<double, unsigned> > rhoenc;
        for (unsigned j = 0; j != nbuckets; ++j) {
          rhoenc.push_back(make_pair(mbucketsi[j]/Dimension::pownu(sqrt((xbucketsi[j] - xi).magnitude2() + softLength2)), j));
        }
        const unsigned ntop = max(1U, unsigned(0.05*nbuckets));
        std::partial_sort(rhoenc.begin(), rhoenc.begin() + ntop, rhoenc.end(), 
                          ComparePairsByFirstElementInDecreasingOrder<pair<double, unsigned> >());

        // For the ntop maxima, sum the contributions from the other centers.  We keep the
        // maximum one.
        for (unsigned k = 0; k != ntop; ++k) {
          const unsigned jmax = rhoenc[k].second;
          const Vector xpcmax = xbucketsi[jmax] - xi;
          const double rpcmax2 = xpcmax.magnitude2();
          double rholocal = 0.0;
          for (unsigned j = 0; j != nbuckets; ++j) {
            const Vector xpc = xbucketsi[j] - xi;
            const double rpc2 = xpc.magnitude2();
            if (rpc2 <= 4.0*rpcmax2 and
                0.75*sqrt(rpcmax2*rpc2) <= xpcmax.dot(xpc)) rholocal += rhoenc[j].first;
          }

          // // Add this nodes local bucket.
          // const LevelKey levelKey = mHomeBuckets(nodeListi, i).first;
          // const CellKey cellKey = mHomeBuckets(nodeListi, i).second;
          // CHECK(levelKey <= mTree.size());
          // typename TreeLevel::const_iterator itr = mTree[levelKey].find(cellKey);
          // CHECK(itr != mTree[levelKey].end());
          // const Cell& localCell = itr->second; // mTree[levelKey][cellKey];
          // rholocal += localCell.Mglobal/Dimension::pownu(sqrt((localCell.xcm - xi).magnitude2() + softLength2));

          // Now check for the global max.
          if (mRhoMax < rholocal) {
            mRhoMax = rholocal;
            mNodeListMax = nodeListi;
            mimax = i;
          }
        }
      }
    }

    // We now have the local maximum effective density for any point.  That's enough
    // to compute the gravitational dynamical time scale.
  }

  // Finalize, and sum up the potential for the extra energy term.
  mExtraEnergy = 0.0;
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mPotential[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {

      // Set the position derivative.
      DxDt(nodeListi, i) = velocity(nodeListi, i);

      // Multiply by G.
      DvDt(nodeListi, i) *= mG;
      mPotential(nodeListi, i) *= mG;

      // Accumluate the package energy as the total gravitational potential.
      mExtraEnergy += 0.5*mass(nodeListi, i)*mPotential(nodeListi, i);
    }
  }

#ifdef ENABLE_MPI
  mExtraEnergy = allReduce(mExtraEnergy, SPHERAL_OP_SUM);

  // Wait until all our sends are complete.
  if (not sendRequests.empty()) {
    vector<MPI_Status> sendStatus(sendRequests.size());
    MPI_Waitall(sendRequests.size(), &(*sendRequests.begin()), &(*sendStatus.begin()));
  }
#endif

}

//------------------------------------------------------------------------------
// Do start of the problem type tasks.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
TreeGravity<Dimension>::
initializeProblemStartup(DataBase<Dimension>& db) {

  // Allocate space for the gravitational potential FieldList.
  mPotential = db.newGlobalFieldList(0.0, "gravitational potential");
}

//------------------------------------------------------------------------------
// Do start of the problem type tasks.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
TreeGravity<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& db,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {
  // We need to make a dry run through setting derivatives and such
  // to set our initial vote on the time step.
  vector<Physics<Dimension>*> packages(1, this);
  this->initialize(0.0, 1.0, db, state, derivs);
  this->evaluateDerivatives(0.0, 1.0, db, state, derivs);
}

//------------------------------------------------------------------------------
// Intialize the package before evaluateDerivatives is called.
// For TreeGravity, this is where we build the current tree.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
TreeGravity<Dimension>::
initialize(const Scalar /*time*/,
           const Scalar /*dt*/,
           const DataBase<Dimension>& db,
           State<Dimension>& state,
           StateDerivatives<Dimension>& /*derivs*/) {

  // For now we're not going to be clever about trying to patch an existing tree,
  // but instead we'll build it from scratch every time.
  mTree = Tree();

  // Access to pertinent fields in the database.
  if (db.numInternalNodes() > 0) {
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
    const size_t numNodeLists = mass.numFields();

    // Determine the box size.
    globalBoundingBox(position, mXmin, mXmax, false);
    CHECK(mXmin.x() <= mXmax.x());
    CHECK(mXmin.y() <= mXmax.y());
    CHECK(mXmin.z() <= mXmax.z());
    mBoxLength = (mXmax - mXmin).maxAbsElement();
    CHECK(mBoxLength >= 0.0);

    // Walk all the internal nodes and add them to the tree.
    for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const size_t n = mass[nodeListi]->numInternalElements();
      for (size_t i = 0; i != n; ++i) {
        this->addNodeToTree(mass(nodeListi, i), position(nodeListi, i), velocity(nodeListi, i));
      }
    }

#ifdef ENABLE_MPI

    // We require the true global center of mass for each cell, which means
    // in parallel we need to exchange the local trees and build up this information.
    // We also build the global total mass in each cell at the same time.
    // Get the processor information.
    const unsigned rank = Process::getRank();
    const unsigned numProcs = Process::getTotalNumberOfProcesses();

    // Pack up the local tree.
    vector<char> localBuffer, buffer;
    this->serialize(mTree, localBuffer);

    // Walk each process, and have them send their local tree info to all other
    // domains.
    unsigned bufSize;
    vector<char>::const_iterator bufItr;
    Tree tree;
    for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {

      // Broadcast the send processor's tree.
      buffer = localBuffer;
      bufSize = buffer.size();
      MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, Communicator::communicator());
      buffer.resize(bufSize);
      MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, Communicator::communicator());
      tree = Tree();
      bufItr = buffer.begin();
      this->deserialize(tree, bufItr, buffer.end());
      CHECK(bufItr == buffer.end());

      // Augment any of our local cell calculations for mass and center of mass.
      if (sendProc != rank) {
        const unsigned nlevels = min(mTree.size(), tree.size());
        for (unsigned ilevel = 0; ilevel != nlevels; ++ilevel) {
          for (typename TreeLevel::iterator myItr = mTree[ilevel].begin();
               myItr != mTree[ilevel].end();
               ++myItr) {
            const CellKey key = myItr->first;
            Cell& cell = myItr->second;
            CHECK(cell.key == key);
            const typename TreeLevel::const_iterator otherItr = tree[ilevel].find(key);
            if (otherItr != tree[ilevel].end()) {
              const Cell& otherCell = otherItr->second;
              cell.xcm = (cell.Mglobal*cell.xcm + otherCell.M*otherCell.xcm)/(cell.Mglobal + otherCell.M);
              cell.vcm = (cell.Mglobal*cell.vcm + otherCell.M*otherCell.vcm)/(cell.Mglobal + otherCell.M);
              cell.Mglobal += otherCell.M;
            }
          }
        }
      }
    }

#endif

    // Set the daughter pointers.
    constructDaughterPtrs(mTree);

    // Make a final pass over the cells and fill in the distance between the 
    // center of mass and the geometric center.
    CellKey ckey, ix, iy, iz;
    double cellsize;
    for (unsigned ilevel = 0; ilevel != mTree.size(); ++ilevel) {
      cellsize = mBoxLength/(1U << ilevel);
      //double cellvol = Dimension::pownu(cellsize);
      for (typename TreeLevel::iterator itr = mTree[ilevel].begin();
           itr != mTree[ilevel].end();
           ++itr) {
        ckey = itr->first;
        Cell& cell = itr->second;
        extractCellIndices(ckey, ix, iy, iz);

        // Update the distance between the cell's center of mass and geometric center.
        cell.rcm2cc2 = (cell.xcm - TreeDimensionTraits<Dimension>::cellCenter(mXmin, ix, iy, iz, cellsize)).magnitude2();
        CHECK(cell.rcm2cc2 < Dimension::nDim*cellsize*cellsize);
      }
    }
  }
  return false;
}

//------------------------------------------------------------------------------
// Vote on a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename TreeGravity<Dimension>::TimeStepType
TreeGravity<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& /*derivs*/,
   const Scalar /*currentTime*/) const {

  // A standard N-body approach -- just take the ratio of softening length/acceleration.
  if (mTimeStepChoice == GravityTimeStepType::AccelerationRatio) {
    const double dt = mftimestep * mDtMinAcc;
    std::stringstream reasonStream;
    reasonStream << "TreeGravity: f*sqrt(L/a) = " << dt << std::endl;
    return TimeStepType(dt, reasonStream.str());

  } else {

    // The following is an implementation (perhaps not exact) of the algorithm 
    // described in Zemp et al., MNRAS, 376, 273-286 (2007).
    // The idea is to come up with the most restrictive dynamical time per particle
    // by find the dominant potentials for each point.
    // Note the computation of rhomax and such is done in evaluateDerivatives.
    CHECK(mRhoMax > 0.0);
    CHECK(mNodeListMax >= 0);
    CHECK(mimax >= 0);
    const double dtDyn = sqrt(1.0/(mG*mRhoMax));
    const double dt = mftimestep * dtDyn;
    const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    std::stringstream reasonStream;
    reasonStream << "TreeGravity: sqrt(1/(G rho)) = sqrt(1/("
                 << mG << " * " << mRhoMax
                 << ")) = " << dt 
                 << " selected for node " << mimax 
                 << " in NodeList " << position[mNodeListMax]->nodeList().name()
                 << std::endl;
    return TimeStepType(dt, reasonStream.str());
  }
}

//------------------------------------------------------------------------------
// extraEnergy
//------------------------------------------------------------------------------
template<typename Dimension>
typename TreeGravity<Dimension>::Scalar 
TreeGravity<Dimension>::
extraEnergy() const {
  return mExtraEnergy;
}

//------------------------------------------------------------------------------
// potential
//------------------------------------------------------------------------------
template<typename Dimension>
const FieldList<Dimension, typename TreeGravity<Dimension>::Scalar>&
TreeGravity<Dimension>::
potential() const {
  return mPotential;
}

//------------------------------------------------------------------------------
// dumpTree
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
TreeGravity<Dimension>::
dumpTree(const bool globalTree) const {
  std::stringstream ss;
  CellKey ix, iy, iz;
#ifdef ENABLE_MPI
  const unsigned numProcs = Process::getTotalNumberOfProcesses();
  const unsigned rank = Process::getRank();
#endif
  unsigned nlevels = mTree.size();
  if (globalTree) nlevels = allReduce(nlevels, SPHERAL_OP_MAX);

  ss << "Tree : nlevels = " << nlevels << "\n";
  for (unsigned ilevel = 0; ilevel != nlevels; ++ilevel) {

    // Gather up the level cells and sort them.
    vector<Cell> cells;
    vector<char> localBuffer;
    cells.reserve(mTree[ilevel].size());
    if (ilevel < mTree.size()) {
      for (typename TreeLevel::const_iterator itr = mTree[ilevel].begin();
           itr != mTree[ilevel].end();
           ++itr) {
        cells.push_back(itr->second);
        this->serialize(itr->second, localBuffer);
      }
    }
#ifdef ENABLE_MPI
    if (globalTree) {
      for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {
        unsigned bufSize = localBuffer.size();
        MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, Communicator::communicator());
        if (bufSize > 0) {
          vector<char> buffer = localBuffer;
          MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, Communicator::communicator());
          if (rank != sendProc) {
            vector<char>::const_iterator itr = buffer.begin();
            while (itr < buffer.end()) {
              cells.push_back(Cell());
              this->deserialize(cells.back(), itr, buffer.end());
            }
          }
        }
      }
    }
#endif
    sort(cells.begin(), cells.end());
    cells.erase(unique(cells.begin(), cells.end()), cells.end());

    ss << "--------------------------------------------------------------------------------\n" 
       << " Level " << ilevel << " : numCells = " << cells.size() << "\n";
    for (typename vector<Cell>::const_iterator itr = cells.begin();
         itr != cells.end();
         ++itr) {
      const Cell& cell = *itr;
      extractCellIndices(cell.key, ix, iy, iz);
      ss << "    Cell key=" << cell.key << " : (ix,iy,iz)=(" << ix << " " << iy << " " << iz << ")\n"
         << "         xcm=" << cell.xcm << " vcm=" << cell.vcm << " rcm2cc=" << sqrt(cell.rcm2cc2) << " M=" << cell.M  << " Mglobal=" << cell.Mglobal << "\n"
         << "         daughters = ( ";
      for (vector<CellKey>::const_iterator ditr = cell.daughters.begin();
           ditr != cell.daughters.end();
           ++ditr) ss << *ditr << " ";
      ss << ")\n"
         << "         nodes = [";
      for (unsigned k = 0; k != cell.masses.size(); ++k) ss << " ("
                                                            << cell.masses[k] << " "
                                                            << cell.positions[k] << " "
                                                            << cell.velocities[k] << ")";
      ss <<" ]\n";
    }
  }
  return ss.str();
}

//------------------------------------------------------------------------------
// dumpTreeStatistics
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
TreeGravity<Dimension>::
dumpTreeStatistics(const bool globalTree) const {
  std::stringstream ss;
#ifdef ENABLE_MPI
  const unsigned numProcs = Process::getTotalNumberOfProcesses();
  const unsigned rank = Process::getRank();
#endif
  unsigned nlevels = mTree.size();
  if (globalTree) nlevels = allReduce(nlevels, SPHERAL_OP_MAX);

  ss << "Tree : nlevels = " << nlevels << "\n";
  for (unsigned ilevel = 0; ilevel != nlevels; ++ilevel) {

    // Gather up the level cells and sort them.
    vector<Cell> cells;
    vector<char> localBuffer;
    cells.reserve(mTree[ilevel].size());
    if (ilevel < mTree.size()) {
      for (typename TreeLevel::const_iterator itr = mTree[ilevel].begin();
           itr != mTree[ilevel].end();
           ++itr) {
        cells.push_back(itr->second);
        this->serialize(itr->second, localBuffer);
      }
    }
#ifdef ENABLE_MPI
    if (globalTree) {
      for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {
        unsigned bufSize = localBuffer.size();
        MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, Communicator::communicator());
        if (bufSize > 0) {
          vector<char> buffer = localBuffer;
          MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, Communicator::communicator());
          if (rank != sendProc) {
            vector<char>::const_iterator itr = buffer.begin();
            while (itr < buffer.end()) {
              cells.push_back(Cell());
              this->deserialize(cells.back(), itr, buffer.end());
            }
          }
        }
      }
    }
#endif
    sort(cells.begin(), cells.end());
    cells.erase(unique(cells.begin(), cells.end()), cells.end());

    // Count up how much of everything we have.
    ss << "--------------------------------------------------------------------------------\n" 
       << " Level " << ilevel << " : numCells = " << cells.size() << "\n";
    unsigned nparticles = 0;
    for (typename vector<Cell>::const_iterator itr = cells.begin();
         itr != cells.end();
         ++itr) {
      const Cell& cell = *itr;
      nparticles += cell.masses.size();
    }
    ss << "         : nparts = " << nparticles << "\n";
  }
  return ss.str();
}

//------------------------------------------------------------------------------
// G
//------------------------------------------------------------------------------
template<typename Dimension>
double
TreeGravity<Dimension>::
G() const {
  return mG;
}

//------------------------------------------------------------------------------
// opening
//------------------------------------------------------------------------------
template<typename Dimension>
double
TreeGravity<Dimension>::
opening() const {
  return sqrt(mOpening2);
}

template<typename Dimension>
void
TreeGravity<Dimension>::
opening(double x) {
  VERIFY(x > 0.0);
  mOpening2 = x*x;
}

//------------------------------------------------------------------------------
// softeningLength
//------------------------------------------------------------------------------
template<typename Dimension>
double
TreeGravity<Dimension>::
softeningLength() const {
  return mSofteningLength;
}

template<typename Dimension>
void
TreeGravity<Dimension>::
softeningLength(double x) {
  VERIFY(x > 0.0);
  mSofteningLength = x;
}

//------------------------------------------------------------------------------
// ftimestep
//------------------------------------------------------------------------------
template<typename Dimension>
double
TreeGravity<Dimension>::
ftimestep() const {
  return mftimestep;
}

template<typename Dimension>
void
TreeGravity<Dimension>::
ftimestep(double x) {
  VERIFY(x > 0.0);
  mftimestep = x;
}

//------------------------------------------------------------------------------
// timeStepChoice
//------------------------------------------------------------------------------
template<typename Dimension>
GravityTimeStepType
TreeGravity<Dimension>::
timeStepChoice() const {
  return mTimeStepChoice;
}

template<typename Dimension>
void
TreeGravity<Dimension>::
timeStepChoice(GravityTimeStepType x) {
  mTimeStepChoice = x;
}

//------------------------------------------------------------------------------
// xmin
//------------------------------------------------------------------------------
template<typename Dimension>
typename TreeGravity<Dimension>::Vector
TreeGravity<Dimension>::
xmin() const {
  return mXmin;
}

//------------------------------------------------------------------------------
// xmax
//------------------------------------------------------------------------------
template<typename Dimension>
typename TreeGravity<Dimension>::Vector
TreeGravity<Dimension>::
xmax() const {
  return mXmax;
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeGravity<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mPotential, pathName + "/potential");
  file.write(mDtMinAcc, pathName + "/pairWiseDtMin");
  file.write(mNodeListMax, pathName + "/nodeListMax");
  file.write(mimax, pathName + "/imax");
  file.write(mRhoMax, pathName + "/rhoMax");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeGravity<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mPotential, pathName + "/potential");
  file.read(mDtMinAcc, pathName + "/pairWiseDtMin");
  file.read(mNodeListMax, pathName + "/nodeListMax");
  file.read(mimax, pathName + "/imax");
  file.read(mRhoMax, pathName + "/rhoMax");
}

//------------------------------------------------------------------------------
// Apply the forces from a tree.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TreeGravity<Dimension>::
applyTreeForces(const Tree& tree,
                const FieldList<Dimension, Scalar>& mass,
                const FieldList<Dimension, Vector>& position,
                FieldList<Dimension, Vector>& DvDt,
                FieldList<Dimension, Scalar>& potential,
                FieldList<Dimension, vector<Scalar> >& interactionMasses,
                FieldList<Dimension, vector<Vector> >& interactionPositions,
                FieldList<Dimension, pair<LevelKey, CellKey> >& homeBuckets,
                TreeGravity<Dimension>::CompletedCellSet& cellsCompleted) const {

  const unsigned numNodeLists = mass.numFields();
  //const unsigned numNodes = mass.numInternalNodes();
  const double boxLength2 = mBoxLength*mBoxLength;
  const double softLength2 = mSofteningLength*mSofteningLength;

  // Prepare the result for the shortest timestep.
  double result = std::numeric_limits<double>::max();

  // Declare variables we're going to need in the loop once.  May help with optimization?
  unsigned nodeListi, i, j, k, ilevel, nremaining;
  double mj, cellsize2, rji2;
  Vector xji, nhat;
  vector<Cell*> remainingCells, newDaughters;
  NodeID inode;

  // We'll always be starting with the daughters of the root level.
  const unsigned numLevels = tree.size();
  if (numLevels >= 1) {
    const Cell& rootCell = tree[0].begin()->second;

    // Walk each internal node.
    typename TreeLevel::const_iterator cellItr;
    for (nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      for (i = 0; i != mass[nodeListi]->numInternalElements(); ++i) {

        // State of node i.
        //mi = mass(nodeListi, i);
        const Vector& xi = position(nodeListi, i);
        Vector& DvDti = DvDt(nodeListi, i);
        Scalar& phii = potential(nodeListi, i);
        inode = NodeID(nodeListi, i);

        // Walk the tree.
        ilevel = 0;
        remainingCells = rootCell.daughterPtrs;
        while ((not remainingCells.empty()) and ++ilevel < numLevels) {
          nremaining = remainingCells.size();
          newDaughters = vector<Cell*>();
          newDaughters.reserve(8*nremaining);
          cellsize2 = boxLength2/(1U << (2*ilevel));

          // Walk each of the current set of Cells.
          for (k = 0; k != nremaining; ++k) {
            const Cell& cell = *remainingCells[k];

            // Can we ignore this cells daughters?
            if (cellsCompleted[inode][ilevel].find(cell.key) == cellsCompleted[inode][ilevel].end()) {
              xji = cell.xcm - xi;
              rji2 = xji.magnitude2();

              // We use Barnes (1994) modified criterion, except for the extra square factor here for efficiency.
              if (rji2 > cellsize2/mOpening2 + cell.rcm2cc2) {      

                // Yep, treat this cells and all of its daughters as a single point.
                nhat = xji.unitVector();
                rji2 += softLength2;
                CHECK(rji2 > 0.0);

                // Increment the acceleration and potential.
                DvDti += cell.Mglobal*TreeDimensionTraits<Dimension>::forceLaw(rji2) * nhat;  // Multiply by G later
                phii -= cell.Mglobal*TreeDimensionTraits<Dimension>::potentialLaw(rji2);   // Multiply by G later
                cellsCompleted[inode][ilevel].insert(cell.key);

                // Add to the set of "interaction" buckets for this point, which are used for
                // choosing the timestep.
                interactionMasses(nodeListi, i).push_back(cell.Mglobal);
                interactionPositions(nodeListi, i).push_back(cell.xcm);

              } else if (cell.daughterPtrs.size() == 0) {

                // This cell represents a leaf (termination of descent.  We just directly
                // add up the node properties of any nodes in the cell.
                CHECK(cell.masses.size() > 0 and
                      cell.positions.size() == cell.masses.size() and
                      cell.velocities.size() == cell.masses.size());
                for (j = 0; j != cell.masses.size(); ++j) {
                  mj = cell.masses[j];
                  const Vector& xj = cell.positions[j];
                  xji = xj - xi;
                  rji2 = xji.magnitude2();

                  if (rji2/softLength2 > 1.0e-10) {           // Screen out self-interaction.
                    nhat = xji.unitVector();
                    rji2 += softLength2;
                    CHECK(rji2 > 0.0);

                    // Increment the acceleration and potential.
                    DvDti += mj*TreeDimensionTraits<Dimension>::forceLaw(rji2) * nhat;  // Multiply by G later
                    phii -= mj*TreeDimensionTraits<Dimension>::potentialLaw(rji2);   // Multiply by G later
                  } else {
                    homeBuckets(nodeListi, i) = make_pair(ilevel, cell.key);
                  }
                }

                // Add this cells contribution to the interaction buckets.
                interactionMasses(nodeListi, i).push_back(cell.Mglobal);
                interactionPositions(nodeListi, i).push_back(cell.xcm);

              } else {

                // We need to walk further down the tree.  Add this cells daughters
                // to the next set.
                copy(cell.daughterPtrs.begin(), cell.daughterPtrs.end(), back_inserter(newDaughters));

              }
            }
          }

          // Update the set of cells to check on the next pass.
          remainingCells = newDaughters;
        }

        // Update the total acceleration based constraint for the time-step.
        const double amag = DvDti.magnitude();
        if (amag > 0.0) result = min(result, sqrt(mSofteningLength / amag));       // Divide by G later
      }
    }
  }
  return result/sqrt(mG);
}

//------------------------------------------------------------------------------
// Serialize a tree to a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeGravity<Dimension>::
serialize(const TreeGravity<Dimension>::Tree& tree,
          std::vector<char>& buffer) const {
  const unsigned nlevels = tree.size();
  packElement(nlevels, buffer);
  for (unsigned ilevel = 0; ilevel != nlevels; ++ilevel) {
    const unsigned ncells = tree[ilevel].size();
    packElement(ncells, buffer);
    for (typename TreeLevel::const_iterator itr = tree[ilevel].begin();
         itr != tree[ilevel].end();
         ++itr) {
      packElement(itr->first, buffer);
      serialize(itr->second, buffer);
    }
  }
}

//------------------------------------------------------------------------------
// Serialize a cell to a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeGravity<Dimension>::
serialize(const TreeGravity<Dimension>::Cell& cell,
          std::vector<char>& buffer) const {
  packElement(cell.M, buffer);
  packElement(cell.Mglobal, buffer);
  packElement(cell.xcm, buffer);
  packElement(cell.vcm, buffer);
  packElement(cell.rcm2cc2, buffer);
  packElement(cell.key, buffer);
  packElement(cell.daughters, buffer);
  packElement(cell.masses, buffer);
  packElement(cell.positions, buffer);
  packElement(cell.velocities, buffer);
}

//------------------------------------------------------------------------------
// Deserialize a tree from a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeGravity<Dimension>::
deserialize(TreeGravity<Dimension>::Tree& tree,
            vector<char>::const_iterator& bufItr,
            const vector<char>::const_iterator& endItr) const {
  unsigned nlevels;
  CellKey key;
  Cell cell;
  unpackElement(nlevels, bufItr, endItr);
  tree.resize(nlevels);
  for (unsigned ilevel = 0; ilevel != nlevels; ++ilevel) {
    unsigned ncells;
    unpackElement(ncells, bufItr, endItr);
    for (unsigned i = 0; i != ncells; ++i) {
      cell.daughters = vector<CellKey>();
      cell.masses = vector<double>();
      cell.positions = vector<Vector>();
      cell.velocities = vector<Vector>();
      unpackElement(key, bufItr, endItr);
      deserialize(cell, bufItr, endItr);
      tree[ilevel][key] = cell;
    }
  }
  constructDaughterPtrs(tree);
}

//------------------------------------------------------------------------------
// Deserialize a cell from a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TreeGravity<Dimension>::
deserialize(TreeGravity<Dimension>::Cell& cell,
            vector<char>::const_iterator& bufItr,
            const vector<char>::const_iterator& endItr) const {
  unpackElement(cell.M, bufItr, endItr);
  unpackElement(cell.Mglobal, bufItr, endItr);
  unpackElement(cell.xcm, bufItr, endItr);
  unpackElement(cell.vcm, bufItr, endItr);
  unpackElement(cell.rcm2cc2, bufItr, endItr);
  unpackElement(cell.key, bufItr, endItr);
  unpackElement(cell.daughters, bufItr, endItr);
  unpackElement(cell.masses, bufItr, endItr);
  unpackElement(cell.positions, bufItr, endItr);
  unpackElement(cell.velocities, bufItr, endItr);
}

}
