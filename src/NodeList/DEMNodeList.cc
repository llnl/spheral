//---------------------------------Spheral++----------------------------------//
// DEMNodeList -- An abstract base class for the NodeLists to represent
//                  fluids.
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "DEM/DEMFieldNames.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DEMNodeList.hh"

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
DEMNodeList<Dimension>::
DEMNodeList(string name,
            const size_t numInternal,
            const size_t numGhost,
            const Scalar hmin,
            const Scalar hmax,
            const Scalar hminratio,
            const Scalar nPerh,
            const Scalar neighborSearchBuffer,
            const size_t maxNumNeighbors):
  NodeList<Dimension>(name, numInternal, numGhost, hmin, hmax, hminratio, nPerh, maxNumNeighbors),
  mNeighborSearchBuffer(neighborSearchBuffer),
  mParticleRadius(DEMFieldNames::particleRadius, *this),
  mCompositeParticleIndex(DEMFieldNames::compositeParticleIndex, *this),
  mUniqueIndex(DEMFieldNames::uniqueIndices, *this){
}

//------------------------------------------------------------------------------
// Set the particle radii
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMNodeList<Dimension>::
particleRadius(const Field<Dimension, typename Dimension::Scalar>& radii) {
  mParticleRadius = radii;
  mParticleRadius.name(DEMFieldNames::particleRadius);
}

//------------------------------------------------------------------------------
// composite indes for conglomerates
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMNodeList<Dimension>::
compositeParticleIndex(const Field<Dimension, int>& ids) {
  mCompositeParticleIndex = ids;
  mCompositeParticleIndex.name(DEMFieldNames::compositeParticleIndex);
}

//------------------------------------------------------------------------------
// Set the unique indices
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMNodeList<Dimension>::
uniqueIndex(const Field<Dimension, size_t>& ids) {
  mUniqueIndex = ids;
  mUniqueIndex.name(DEMFieldNames::uniqueIndices);
}
//------------------------------------------------------------------------------
// Calculate and return the volume per node.
//------------------------------------------------------------------------------
//template<typename Dimension>
//void
//DEMNodeList<Dimension>::
//volume(Field<Dimension, typename Dimension::Scalar>& field) const {
//  REQUIRE(field.nodeListPtr() == this);
//  const Field<Dimension, Scalar>& nodeMass = this->mass();
//  const Field<Dimension, Scalar>& nodeDensity = this->angularVelocity();
//  for (size_t i = 0; i != this->numNodes(); ++i) field[i] = nodeMass[i]*safeInv(nodeDensity[i]);
//  field.name(HydroFieldNames::volume);
//}

//------------------------------------------------------------------------------
// Calculate and return the linear momentum per node.
//------------------------------------------------------------------------------
// template<typename Dimension>
// void
// DEMNodeList<Dimension>::
// linearMomentum(Field<Dimension, typename Dimension::Vector>& field) const {
//   REQUIRE(field.nodeListPtr() == this);
//   const Field<Dimension, Scalar>& nodeMass = this->mass();
//   const Field<Dimension, Vector>& nodeVelocity = this->velocity();
//   for (size_t i = 0; i != this->numNodes(); ++i) field[i] = nodeMass[i]*nodeVelocity[i];
//   field.name(HydroFieldNames::linearMomentum);
// }

//------------------------------------------------------------------------------
// Calculate and return the total energy per node.
//------------------------------------------------------------------------------
// template<typename Dimension>
// void
// DEMNodeList<Dimension>::
// totalEnergy(Field<Dimension, typename Dimension::Scalar>& field) const {
//   REQUIRE(field.nodeListPtr() == this);
//   const Field<Dimension, Scalar>& nodeMass = this->mass();
//   const Field<Dimension, Vector>& nodeVelocity = this->velocity();
//   const Field<Dimension, Scalar>& nodeEps = this->specificThermalEnergy();
//   for (size_t i = 0; i != this->numNodes(); ++i) 
//     field[i] = nodeMass[i]*(0.5*nodeVelocity[i].magnitude2() + nodeEps[i]);
//   field.name(HydroFieldNames::totalEnergy);
// }


//------------------------------------------------------------------------------
// when problem starts set our equilibrium overlap
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMNodeList<Dimension>::
setHfieldFromParticleRadius(const size_t startUniqueIndex){
  
  const auto kernelExtent = this->neighbor().kernelExtent();
  const auto& radius = this->particleRadius();
  const auto& uniqId = this->uniqueIndex();
  auto& Hfield = this->Hfield();
  const auto ni = this->numInternalNodes();

#pragma omp parallel for
  for (auto i = 0u; i < ni; ++i) {
    const auto ui = uniqId[i];
    if(ui >= startUniqueIndex) {
      const auto Ri = radius[i];
      auto& Hi = Hfield[i];
      const auto hInv = safeInv(2.0 * Ri * (1.0+mNeighborSearchBuffer)/kernelExtent);
      Hi = SymTensor::one() * hInv;
    }
  }
}


//------------------------------------------------------------------------------
// Dump the current state of the NodeList to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMNodeList<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  // Dump the ancestor class.
  NodeList<Dimension>::dumpState(file, pathName);

  // Dump each of the internal fields of the DEMNodeList.
  file.write(mParticleRadius, pathName + "/particleRadius");
  file.write(mCompositeParticleIndex, pathName + "/compositeParticleIndex");
  file.write(mUniqueIndex, pathName + "/uniqueIndex");
}

//------------------------------------------------------------------------------
// Restore the state of the NodeList from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMNodeList<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  // Restore the ancestor class.
  NodeList<Dimension>::restoreState(file, pathName);

  // Restore each of the internal fields of the DEMNodeList.
  file.read(mParticleRadius, pathName + "/particleRadius");
  file.read(mCompositeParticleIndex, pathName + "/compositeParticleIndex");
  file.read(mUniqueIndex, pathName + "/uniqueIndex");
}  

}
