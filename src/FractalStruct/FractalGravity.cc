//---------------------------------Spheral++----------------------------------//
// FractalGravity implementation.
//
//! \author $Author: mikeowen $
//! \version $Revision: 4012 $
//! \date $Date: 2011-02-20 10:58:32 -0800 (Sun, 20 Feb 2011) $
//----------------------------------------------------------------------------//

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "FractalGravity.hh"
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
#include "Distributed/Communicator.hh"

#include "libs.hh"
#include "classes.hh"
#include "headers.hh"

#include <mpi.h>

namespace Spheral {
namespace GravitySpace {

using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
FractalGravity::
FractalGravity(const double G,
               const Vector& xmin,
               const Vector& xmax,
               const bool periodic,
               const unsigned gridLength,
               const double maxDeltaVelocity):
  mG(G),
  mXmin(xmin),
  mXmax(xmax),
  mPeriodic(periodic),
  mGridLength(gridLength),
  mMaxDeltaVelocityFactor(maxDeltaVelocity),
  mPotential(FieldSpace::FieldStorageType::CopyFields),
  mExtraEnergy(0.0),
  mOldMaxAcceleration(0.0),
  mOldMaxVelocity(0.0),
  mFractalMemoryPtr(NULL),
  mFractalComm() {

  // Decide how many MPI ranks we're going to use, and create the communicator.
  // For now we assume we're modeling a cube, so just make the cpu ranks along each direction equal.
  const auto mxFractalNodes = std::max(1, int(pow(double(Process::getTotalNumberOfProcesses()), 1.0/3.0)));
  const auto maxProc = mxFractalNodes*mxFractalNodes*mxFractalNodes;
  CHECK(maxProc <= Process::getTotalNumberOfProcesses());
  MPI_Comm_split(Communicator::communicator(),
                 (Process::getRank() < maxProc ? 1 : MPI_UNDEFINED),
                 Process::getRank(),
                 &mFractalComm);

  // Call Jens' setup method.
  if (mPeriodic) {
    mFractalMemoryPtr = FractalSpace::FractalGravityFirstTime(mPeriodic,
                                                              mFractalComm,
                                                              mGridLength,
                                                              mxFractalNodes,
                                                              mxFractalNodes,
                                                              mxFractalNodes,
                                                              "",
                                                              "blago");
  } else {
    mFractalMemoryPtr = FractalSpace::FractalGravityIsolatedFirstTime(mFractalComm,
                                                                      mGridLength,
                                                                      mxFractalNodes,
                                                                      mxFractalNodes,
                                                                      mxFractalNodes,
                                                                      "",
                                                                      "blago");
  }
  FractalSpace::fractal_memory_setup(mFractalMemoryPtr);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
FractalGravity::
~FractalGravity() {
  FractalSpace::fractal_memory_content_delete(mFractalMemoryPtr);
  FractalSpace::fractal_memory_delete(mFractalMemoryPtr);
}

//------------------------------------------------------------------------------
// Evaluate the forces and return our time derivatives.
//------------------------------------------------------------------------------
void 
FractalGravity::
evaluateDerivatives(const Dim<3>::Scalar time,
                    const Dim<3>::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dim<3> >& state,
                    StateDerivatives<Dim<3> >& derivs) const {
  using namespace NodeSpace;

  // Access to pertinent fields in the database.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto numNodeLists = mass.numFields();
  const auto numNodes = mass.numInternalNodes();

  // Get the acceleration and position change vectors we'll be modifying.
  auto DxDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  auto DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);

  // Make Fractal friendly versions of the fields.
  vector<double> m(numNodes), xpos(numNodes), ypos(numNodes), zpos(numNodes);
  {
    auto j = 0;
    for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
      const auto n = mass[nodeListi]->numInternalElements();
      for (auto i = 0; i < n; ++i) {
        m[j] = mass(nodeListi, i);
        xpos[j] = position(nodeListi, i).x();
        ypos[j] = position(nodeListi, i).y();
        zpos[j] = position(nodeListi, i).z();
        ++j;
      }
    }
    CHECK(j == numNodes);
  }

  // Zero out the total gravitational potential energy.
  mExtraEnergy = 0.0;

  // Call Fractal's innards to do the force evaluations.
  mFractalMemoryPtr->setNumberParticles(numNodes);
  FractalSpace::fractal_create(mFractalMemoryPtr);
  FractalSpace::add_particles(mFractalMemoryPtr,
                              0,                     // first particle index
                              numNodes,              // last particle index
                              xpos,                  // the coordinates
                              ypos,
                              zpos,
                              m);                    // masses
  vector<double> fxmin(3), fxmax(3),
                 fxmin_in(mXmin.begin(), mXmin.end()),
                 fxmax_in(mXmax.begin(), mXmax.end());
  FractalSpace::FractalCube(mFractalMemoryPtr,
                            (mPeriodic ? 0.0 : 1.0), // SHRINK: 0=>use supplied bounds, 1=>compute bounds
                            fxmin_in,                // Input box min coordinate
                            fxmax_in,                // Input box max coordinate
                            fxmin,                   // Output box min coordinate
                            fxmax);                  // Output box max coordinate
  FractalSpace::balance_by_particles(mFractalMemoryPtr, true);
  FractalSpace::DoFractalGravity(mFractalMemoryPtr);

  // Extract the potential and accelerations.
  vector<double> phi(numNodes), accx(numNodes), accy(numNodes), accz(numNodes);
  FractalSpace::get_field(mFractalMemoryPtr,
                          0,                         // first particle index
                          numNodes,                  // last particle index
                          mG,                        // gravitational constant
                          fxmin,                     // min box coordinates
                          fxmax,                     // max box coordinates
                          phi,                       // potential
                          accx,                      // acceleration components
                          accy,
                          accz);
  FractalSpace::fractal_delete(mFractalMemoryPtr);

  // Fill Spheral's fields with the computed values.
  {
    auto j = 0;
    for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
      const auto n = mass[nodeListi]->numInternalElements();
      for (auto i = 0; i < n; ++i) {
        DvDt(nodeListi, i) += Vector(accx[j], accy[j], accz[j]);
        mPotential(nodeListi, i) = phi[j];
        mExtraEnergy += 0.5*mass(nodeListi, i)*phi[j];
        ++j;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Do start of the problem type tasks.
//------------------------------------------------------------------------------
void 
FractalGravity::
initializeProblemStartup(DataBase<Dimension>& db) {

  // Allocate space for the gravitational potential FieldList.
  mPotential = db.newGlobalFieldList(0.0, "gravitational potential");

}

//------------------------------------------------------------------------------
// Vote on a time step.  We should fill in a sqrt(G/rho) type thing here!
//------------------------------------------------------------------------------
FractalGravity::TimeStepType
FractalGravity::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {

  // The maximum change in our velocity for the next time cycle is given 
  // by the mMaxDeltaVelocityFactor plus the max velocity and the max 
  // accelation from the last cycle.
  const double deltat = (mOldMaxVelocity*mOldMaxAcceleration/(mOldMaxAcceleration*mOldMaxAcceleration + 1.0e-10)) * mMaxDeltaVelocityFactor;

  stringstream reasonStream;
  reasonStream << "velocity: " << mOldMaxVelocity
               << ", acceleration: " << mOldMaxAcceleration
               << "dt = f*v/a: " << deltat
               << ends;
  return TimeStepType(deltat, reasonStream.str());
}

//------------------------------------------------------------------------------
// extraEnergy
//------------------------------------------------------------------------------
FractalGravity::Scalar 
FractalGravity::
extraEnergy() const {
  return mExtraEnergy;
}

//------------------------------------------------------------------------------
// potential
//------------------------------------------------------------------------------
const FieldList<Dim<3>, FractalGravity::Scalar>&
FractalGravity::
potential() const {
  return mPotential;
}

//------------------------------------------------------------------------------
// G
//------------------------------------------------------------------------------
double
FractalGravity::
G() const {
  return mG;
}

//------------------------------------------------------------------------------
// xmin
//------------------------------------------------------------------------------
FractalGravity::Vector
FractalGravity::
xmin() const {
  return mXmin;
}

//------------------------------------------------------------------------------
// xmax
//------------------------------------------------------------------------------
FractalGravity::Vector
FractalGravity::
xmax() const {
  return mXmax;
}

//------------------------------------------------------------------------------
// periodic
//------------------------------------------------------------------------------
bool
FractalGravity::
periodic() const {
  return mPeriodic;
}

//------------------------------------------------------------------------------
// gridLength
//------------------------------------------------------------------------------
unsigned
FractalGravity::
gridLength() const {
  return mGridLength;
}

}
}
