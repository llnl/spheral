//---------------------------------Spheral++----------------------------------//
// NBodyGravity implementation.
//
//! \author $Author: mikeowen $
//! \version $Revision: 4012 $
//! \date $Date: 2011-02-20 10:58:32 -0800 (Sun, 20 Feb 2011) $
//----------------------------------------------------------------------------//

#include "NBodyGravity.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "Utilities/DBC.hh"
#include "Material/PhysicalConstants.hh"
#include "Utilities/packElement.hh"

#include "config.hh" // for ENABLE_MPI

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
template<typename Dimension> struct GravityDimensionTraits;

//..............................................................................
// 1D
//..............................................................................
template<>
struct GravityDimensionTraits<Dim<1>> {
  static double forceLaw(const double) { return 1.0; }
  static double potentialLaw(const double) { return 1.0; }  // Not sure about this one
};

//..............................................................................
// 2D
//..............................................................................
template<>
struct GravityDimensionTraits<Dim<2>> {
  static double forceLaw(const double r2) { return 1.0/sqrt(r2); }
  static double potentialLaw(const double r2) { return log(sqrt(r2)); }
};

//..............................................................................
// 3D
//..............................................................................
template<>
struct GravityDimensionTraits<Dim<3>> {
  static double forceLaw(const double r2) { return 1.0/r2; }
  static double potentialLaw(const double r2) { return 1.0/sqrt(r2); }
};

}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template <typename Dimension>
NBodyGravity<Dimension>::
NBodyGravity(const double plummerSofteningLength,
             const double maxDeltaVelocity,
             const double G,
             const bool compatibleVelocityUpdate):
  mPotential(FieldStorageType::CopyFields),
  mPotential0(FieldStorageType::CopyFields),
  mVel02(FieldStorageType::CopyFields),
  mExtraEnergy(0.0),
  mMaxDeltaVelocityFactor(maxDeltaVelocity),
  mSofteningLength(plummerSofteningLength),
  mG(G),
  mCompatibleVelocityUpdate(compatibleVelocityUpdate) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template <typename Dimension>
NBodyGravity<Dimension>::
~NBodyGravity() {
}

//------------------------------------------------------------------------------
// Register some extra state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NBodyGravity<Dimension>::
registerState(DataBase<Dimension >& dataBase,
              State<Dimension >& state) {

  GenericBodyForce<Dimension >::registerState(dataBase, state);
  state.enroll(mPotential);
}

//------------------------------------------------------------------------------
// Derivatives
//------------------------------------------------------------------------------
template <typename Dimension>
void 
NBodyGravity<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Find the square of the Plummer softening length.
  //auto softeningLength2 = mSofteningLength * mSofteningLength;

  // Access to pertinent fields in the database.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);

  // Get the acceleration and position change vectors we'll be modifying.
  auto DxDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  auto DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);

  // Zero out the total gravitational potential energy.
  mExtraEnergy = 0.0;
  mPotential = 0.0;
  mOldMaxAcceleration = 1.0e-10;
  mOldMaxVelocity = 0.0;

  // Pack up the local particle info.
  vector<char> localBuffer, buffer;
  this->serialize(mass, position, localBuffer);

#ifdef ENABLE_MPI
  // Get the processor information.
  const unsigned rank = Process::getRank();
  const unsigned numProcs = Process::getTotalNumberOfProcesses();

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
#endif

  // Add our local contributions.
  vector<Scalar> otherMass;
  vector<Vector> otherPosition;
  this->deserialize(localBuffer, otherMass, otherPosition);
  this->applyPairForces(otherMass, otherPosition, position, DvDt, mPotential);

#ifdef ENABLE_MPI
  // Now walk the other processes and get their contributions.
  unsigned bufSize;
  MPI_Status recvStatus;
  for (unsigned otherProc = 0; otherProc != numProcs; ++otherProc) {
    if (otherProc != rank) {
      MPI_Recv(&bufSize, 1, MPI_UNSIGNED, otherProc, 1, Communicator::communicator(), &recvStatus);
      if (bufSize > 0) {
        buffer = vector<char>(bufSize);
        MPI_Recv(&buffer.front(), bufSize, MPI_CHAR, otherProc, 2, Communicator::communicator(), &recvStatus);
        this->deserialize(buffer, otherMass, otherPosition);
        this->applyPairForces(otherMass, otherPosition,position, DvDt, mPotential);
      }
    }
  }
#endif

  // Finalize our stuff.
  //const unsigned numNodeLists = dataBase.numNodeLists();
  unsigned ifield = 0;
  for (auto iitr = dataBase.fluidNodeListBegin();
       iitr != dataBase.fluidNodeListEnd();
       ++iitr, ++ifield) {
    const auto& nodeListi = **iitr;
    const auto n = nodeListi.numInternalNodes();
    for (auto i = 0u; i != n; ++i) {

      // Set the position derivative.
      DxDt(ifield, i) = velocity(ifield, i);

      // Multiply by G.
      DvDt(ifield, i) *= mG;
      mPotential(ifield, i) *= mG;

      // Accumluate the package energy as the total gravitational potential.
      // mExtraEnergy += mPotential(ifield, i);
      mExtraEnergy += 0.5*mass(ifield, i)*mPotential(ifield, i);

      // Capture the maximum acceleration and velocity magnitudes.
      //const auto accelMagnitude = DvDt(ifield, i).magnitude();
      mOldMaxAcceleration = std::max(mOldMaxAcceleration, DvDt(ifield, i).magnitude());
      mOldMaxVelocity = std::max(mOldMaxVelocity, velocity(ifield, i).magnitude());
    }
  }

  mExtraEnergy = allReduce(mExtraEnergy, SPHERAL_OP_SUM);
  mOldMaxAcceleration = allReduce(mOldMaxAcceleration, SPHERAL_OP_MAX);
  mOldMaxVelocity = allReduce(mOldMaxVelocity, SPHERAL_OP_MAX);

#ifdef ENABLE_MPI
  // Wait until all our sends are complete.
  vector<MPI_Status> sendStatus(sendRequests.size());
  MPI_Waitall(sendRequests.size(), &(*sendRequests.begin()), &(*sendStatus.begin()));
#endif
}

//------------------------------------------------------------------------------
// Problem startup
//------------------------------------------------------------------------------
template <typename Dimension>
void 
NBodyGravity<Dimension>::
initializeProblemStartup(DataBase<Dimension>& db) {

  // Allocate space for the gravitational potential FieldList.
  mPotential = db.newGlobalFieldList(0.0, "gravitational potential");
  mPotential0 = db.newGlobalFieldList(0.0, "gravitational potential 0");
  mVel02 = db.newGlobalFieldList(0.0, "vel0 square");
  mPotential.copyFields();
  mPotential0.copyFields();
  mVel02.copyFields();
}

//------------------------------------------------------------------------------
// Do start of the problem type tasks.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
NBodyGravity<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& db,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {
  vector<Physics<Dimension>*> packages(1, this);
  this->initialize(0.0, 1.0, db, state, derivs);
  this->evaluateDerivatives(0.0, 1.0, db, state, derivs);
}

//------------------------------------------------------------------------------
// Pre-step initialization
//------------------------------------------------------------------------------
template <typename Dimension>
void 
NBodyGravity<Dimension>::
preStepInitialize(const DataBase<Dimension>& /*dataBase*/, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {

  if (mCompatibleVelocityUpdate) {

    // Copy the starting potential.
    mPotential0 = mPotential;
    mPotential0.copyFields();
  
    // Take a snapshot of the starting velocity^2 (for KE0).
    const auto vel = state.fields(HydroFieldNames::velocity, Vector::zero);
    const auto numNodeLists = vel.numFields();
    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      const auto n = vel[nodeListi]->numInternalElements();
      for (auto i = 0u; i < n; ++i) {
        mVel02(nodeListi, i) = vel(nodeListi, i).magnitude2();
      }
    }
  }
}

//------------------------------------------------------------------------------
// Post-step finalizations.
//------------------------------------------------------------------------------
template <typename Dimension>
void
NBodyGravity<Dimension>::
finalize(const Scalar time, 
         const Scalar dt,
         DataBase<Dimension>& dataBase, 
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {

  // Augment the kinetic energy to exactly balance the potential energy change.
  if (mCompatibleVelocityUpdate) {

    this->evaluateDerivatives(time, dt, dataBase, state, derivs);

    // Assume mPotential holds the correct end-of-step potential at this time.
    // Correct for Verlet integrator.
    const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
    const auto mass = state.fields(HydroFieldNames::mass, 0.0);
    auto       vel = state.fields(HydroFieldNames::velocity, Vector::zero);
    const auto numNodeLists = vel.numFields();
    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      const auto n = vel[nodeListi]->numInternalElements();
      for (auto i = 0u; i < n; ++i) {
        const auto vi1t2 = mVel02(nodeListi, i) - 2.0*(mPotential(nodeListi, i) - mPotential0(nodeListi, i));
        CHECK(vi1t2 >= 0.0);
        const Scalar f = sqrt(vi1t2*safeInvVar(vel(nodeListi, i).magnitude2()));
        // cerr << " --> " << i << " " << f << " "  << (mPotential(nodeListi, i) - mPotential0(nodeListi, i)) << endl;
        // cerr << "     " << mPotential(nodeListi, i) << " " << (mG*mass(nodeListi, 1)/sqrt((pos(nodeListi, 0) - pos(nodeListi, 1)).magnitude2() + mSofteningLength*mSofteningLength)) << endl;
        vel(nodeListi, i) *= f;
      }
    }
  }

}

//------------------------------------------------------------------------------
// Timestep vote
//------------------------------------------------------------------------------
template <typename Dimension>
typename NBodyGravity<Dimension>::TimeStepType
NBodyGravity<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/, 
   const State<Dimension>& /*state*/,
   const StateDerivatives<Dimension>& /*derivs*/,
   const Scalar /*currentTime*/) const {

  // The maximum change in our velocity for the next time cycle is given 
  // by the mMaxDeltaVelocityFactor plus the max velocity and the max 
  // accelation from the last cycle.
  const double deltat = std::min((mOldMaxVelocity*mOldMaxAcceleration/(mOldMaxAcceleration*mOldMaxAcceleration + 1.0e-10)) * mMaxDeltaVelocityFactor,
                                 sqrt(2.0*mMaxDeltaVelocityFactor*mSofteningLength*safeInvVar(mOldMaxAcceleration, 1.0e-10)));

  std::stringstream reasonStream;
  reasonStream << "velocity: " << mOldMaxVelocity
               << ", acceleration: " << mOldMaxAcceleration
               << "dt = f*v/a: " << deltat
               << std::endl;
  return TimeStepType(deltat, reasonStream.str());
}

//------------------------------------------------------------------------------
// Package energy
//------------------------------------------------------------------------------
template <typename Dimension>
typename NBodyGravity<Dimension>::Scalar 
NBodyGravity<Dimension>::
extraEnergy() const {
  return mExtraEnergy;
}

//------------------------------------------------------------------------------
// Access the potential field
//------------------------------------------------------------------------------
template <typename Dimension>
const FieldList<Dimension, typename NBodyGravity<Dimension>::Scalar>&
NBodyGravity<Dimension>::
potential() const {
  return mPotential;
}

//------------------------------------------------------------------------------
// valid
//------------------------------------------------------------------------------
template <typename Dimension>
bool 
NBodyGravity<Dimension>::
valid() const {
  // We're always valid, sir!  (This is crap, but we can make no other 
  // assumptions right now.)
  return true;
}

//------------------------------------------------------------------------------
// G
//------------------------------------------------------------------------------
template <typename Dimension>
double
NBodyGravity<Dimension>::
G() const {
  return mG;
}

//------------------------------------------------------------------------------
// softening length
//------------------------------------------------------------------------------
template <typename Dimension>
double
NBodyGravity<Dimension>::
softeningLength() const {
  return mSofteningLength;
}

template <typename Dimension>
void
NBodyGravity<Dimension>::
softeningLength(const double x) {
  VERIFY(x >= 0.0);
  mSofteningLength = x;
}

//------------------------------------------------------------------------------
// compatible velocity
//------------------------------------------------------------------------------
template <typename Dimension>
bool
NBodyGravity<Dimension>::
compatibleVelocityUpdate() const {
  return mCompatibleVelocityUpdate;
}

template <typename Dimension>
void
NBodyGravity<Dimension>::
compatibleVelocityUpdate(const bool x) {
  mCompatibleVelocityUpdate = x;
}

//------------------------------------------------------------------------------
// Accumulate the pair-forces from the given points.
//------------------------------------------------------------------------------
template <typename Dimension>
void 
NBodyGravity<Dimension>::
applyPairForces(const std::vector<Scalar>& otherMass,
                const std::vector<Vector>& otherPosition,
                const FieldList<Dimension, Vector>& position,
                FieldList<Dimension, Vector>& DvDt,
                FieldList<Dimension, Scalar>& potential) const {


  const unsigned numNodeLists = position.numFields();
  const unsigned nother = otherMass.size();
  CHECK(otherPosition.size() == nother);

  // Find the square of the Plummer softening length.
  auto softeningLength2 = mSofteningLength * mSofteningLength;
  auto rmin = 1e-10*mSofteningLength;

  // Some scratch variables.
  unsigned nodeListi, n, i, j;
  Scalar r2;
  Vector r, rHat;

  // Loop over each particle...
  for (nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    n = position[nodeListi]->numInternalElements();
    for (i = 0; i != n; ++i) {
      auto& DvDti = DvDt(nodeListi, i);
      auto& phii = potential(nodeListi, i);
    
      // Walk the other points.
      for (j = 0; j != nother; ++j) {
        r = position(nodeListi, i) - otherPosition[j];
        r2 = r.magnitude2();

        // Particles can't self-interact, silly!
        if (r2 > rmin) {
          r2 += softeningLength2;
          rHat = r.unitVector();

          // Force
          DvDti -= otherMass[j] * rHat * GravityDimensionTraits<Dimension>::forceLaw(r2);  // Multiply by G later

          // Potential
          phii -= otherMass[j] * GravityDimensionTraits<Dimension>::potentialLaw(r2);      // Multiply by G later
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// Serialize
//------------------------------------------------------------------------------
template <typename Dimension>
void
NBodyGravity<Dimension>::
serialize(const FieldList<Dimension, typename Dimension::Scalar>& mass,
          const FieldList<Dimension, typename Dimension::Vector>& position,
          std::vector<char>& buffer) const {
  const unsigned n = mass.numInternalElements();
  CHECK(position.numInternalElements() == n);
  packElement(n, buffer);
  const unsigned numFields = mass.numFields();
  for (unsigned ifield = 0; ifield != numFields; ++ifield) {
    const unsigned n = mass[ifield]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      packElement(mass(ifield, i), buffer);
      packElement(position(ifield, i), buffer);
    }
  }
}

//------------------------------------------------------------------------------
// Deserialize
//------------------------------------------------------------------------------
template <typename Dimension>
void
NBodyGravity<Dimension>::
deserialize(const std::vector<char>& buffer,
            std::vector<typename Dimension::Scalar>& mass,
            std::vector<typename Dimension::Vector>& position) const {
  auto bufItr = buffer.begin();
  unsigned n;
  unpackElement(n, bufItr, buffer.end());
  mass = vector<Scalar>(n);
  position = vector<Vector>(n);
  for (unsigned i = 0; i != n; ++i) {
    unpackElement(mass[i], bufItr, buffer.end());
    unpackElement(position[i], bufItr, buffer.end());
  }
  CHECK(bufItr == buffer.end());
}

} // end namespace Spheral

