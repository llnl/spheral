//---------------------------------Spheral++----------------------------------//
// FluidNodeList -- An abstract base class for the NodeLists to represent
//                  fluids.
//
// Created by JMO, Sat Sep 18 10:50:42 PDT 1999
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Material/EquationOfState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "Kernel/TableKernel.hh"
#include "Field/FieldList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/safeInv.hh"
#include "FluidNodeList.hh"

using std::vector;
using std::list;
using std::string;
using std::cerr;
using std::endl;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with given EOS object, along with optional numInternal nodes,
// numGhost nodes, and name.
//------------------------------------------------------------------------------
template<typename Dimension>
FluidNodeList<Dimension>::
FluidNodeList(string name,
              EquationOfState<Dimension>& eos,
              const size_t numInternal,
              const size_t numGhost,
              const Scalar hmin,
              const Scalar hmax,
              const Scalar hminratio,
              const Scalar nPerh,
              const size_t maxNumNeighbors,
              const Scalar rhoMin,
              const Scalar rhoMax):
  NodeList<Dimension>(name, numInternal, numGhost, hmin, hmax, hminratio, nPerh, maxNumNeighbors),
  mRhoMin(rhoMin),
  mRhoMax(rhoMax),
  mMassDensity(HydroFieldNames::massDensity, *this),
  mSpecificThermalEnergy(HydroFieldNames::specificThermalEnergy, *this),
  mEosPtr(&eos) {
}

//------------------------------------------------------------------------------
// Set the mass density field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FluidNodeList<Dimension>::
massDensity(const Field<Dimension, typename Dimension::Scalar>& rho) {
  mMassDensity = rho;
  mMassDensity.name(HydroFieldNames::massDensity);
}

//------------------------------------------------------------------------------
// Set the specific thermal energy field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FluidNodeList<Dimension>::
specificThermalEnergy(const Field<Dimension, typename Dimension::Scalar>& eps) {
  mSpecificThermalEnergy = eps;
  mSpecificThermalEnergy.name(HydroFieldNames::specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Compute the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FluidNodeList<Dimension>::
pressure(Field<Dimension, typename Dimension::Scalar>& field) const {
  REQUIRE(field.nodeListPtr() == this);
  mEosPtr->setPressure(field, mMassDensity, mSpecificThermalEnergy);
  field.name(HydroFieldNames::pressure);
}

//------------------------------------------------------------------------------
// Calculate and return the temperature field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FluidNodeList<Dimension>::
temperature(Field<Dimension, typename Dimension::Scalar>& field) const {
  REQUIRE(field.nodeListPtr() == this);
  mEosPtr->setTemperature(field, mMassDensity, mSpecificThermalEnergy);
  field.name(HydroFieldNames::temperature);
}

//------------------------------------------------------------------------------
// Calculate and return the sound speed field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FluidNodeList<Dimension>::
soundSpeed(Field<Dimension, typename Dimension::Scalar>& field) const {
  REQUIRE(field.nodeListPtr() == this);
  mEosPtr->setSoundSpeed(field, mMassDensity, mSpecificThermalEnergy);
  field.name(HydroFieldNames::soundSpeed);
}

//------------------------------------------------------------------------------
// Calculate and return the volume per node.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FluidNodeList<Dimension>::
volume(Field<Dimension, typename Dimension::Scalar>& field) const {
  REQUIRE(field.nodeListPtr() == this);
  const Field<Dimension, Scalar>& nodeMass = this->mass();
  const Field<Dimension, Scalar>& nodeDensity = this->massDensity();
  for (size_t i = 0; i != this->numNodes(); ++i) field[i] = nodeMass[i]*safeInv(nodeDensity[i]);
  field.name(HydroFieldNames::volume);
}

//------------------------------------------------------------------------------
// Calculate and return the linear momentum per node.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FluidNodeList<Dimension>::
linearMomentum(Field<Dimension, typename Dimension::Vector>& field) const {
  REQUIRE(field.nodeListPtr() == this);
  const Field<Dimension, Scalar>& nodeMass = this->mass();
  const Field<Dimension, Vector>& nodeVelocity = this->velocity();
  for (size_t i = 0; i != this->numNodes(); ++i) field[i] = nodeMass[i]*nodeVelocity[i];
  field.name(HydroFieldNames::linearMomentum);
}

//------------------------------------------------------------------------------
// Calculate and return the total energy per node.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FluidNodeList<Dimension>::
totalEnergy(Field<Dimension, typename Dimension::Scalar>& field) const {
  REQUIRE(field.nodeListPtr() == this);
  const Field<Dimension, Scalar>& nodeMass = this->mass();
  const Field<Dimension, Vector>& nodeVelocity = this->velocity();
  const Field<Dimension, Scalar>& nodeEps = this->specificThermalEnergy();
  for (size_t i = 0; i != this->numNodes(); ++i) 
    field[i] = nodeMass[i]*(0.5*nodeVelocity[i].magnitude2() + nodeEps[i]);
  field.name(HydroFieldNames::totalEnergy);
}

//------------------------------------------------------------------------------
// Dump the current state of the NodeList to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FluidNodeList<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  // Dump the ancestor class.
  NodeList<Dimension>::dumpState(file, pathName);

  // Dump each of the internal fields of the FluidNodeList.
  file.write(mMassDensity, pathName + "/massDensity");
  file.write(mSpecificThermalEnergy, pathName + "/specificThermalEnergy");
}

//------------------------------------------------------------------------------
// Restore the state of the NodeList from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FluidNodeList<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  // Restore the ancestor class.
  NodeList<Dimension>::restoreState(file, pathName);

  // Restore each of the internal fields of the FluidNodeList.
  file.read(mMassDensity, pathName + "/massDensity");
  file.read(mSpecificThermalEnergy, pathName + "/specificThermalEnergy");
}  

}
