//---------------------------------Spheral++----------------------------------//
// DataBase -- The central point to store NodeLists and fields, to serve the 
//             Spheral++ physics modules.
//
// Created by JMO, Sun Feb  6 13:44:49 PST 2000
//----------------------------------------------------------------------------//
#include "DataBase/DataBase.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/GeometryRegistrar.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Kernel/TableKernel.hh"
#include "Material/EquationOfState.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/safeInv.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/SpheralMessage.hh"
#include "Distributed/allReduce.hh"
#include "Distributed/Communicator.hh"
#include "Utilities/DBC.hh"

#ifdef SPHERAL_ENABLE_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <memory>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
DataBase<Dimension>::
DataBase():
  DataBaseView<Dimension>(),
  mNodeListPtrs(),
  mFluidNodeListPtrs(),
  mFluidNodeListAsNodeListPtrs(),
  mSolidNodeListPtrs(),
  mSolidNodeListAsNodeListPtrs(),
  mDEMNodeListPtrs(),
  mDEMNodeListAsNodeListPtrs(),
  mConnectivityMapPtr(nullptr),
  mNodeListViewStorage(),
  mFluidNodeListViewStorage(),
  mSolidNodeListViewStorage(),
  mDEMNodeListViewStorage() {
}

//------------------------------------------------------------------------------
// Global numbers of nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
size_t
DataBase<Dimension>::globalNumInternalNodes() const {
  return allReduce(this->numInternalNodes(), SPHERAL_OP_SUM);
}

template<typename Dimension>
size_t
DataBase<Dimension>::globalNumGhostNodes() const {
  return allReduce(this->numGhostNodes(), SPHERAL_OP_SUM);
}

template<typename Dimension>
size_t
DataBase<Dimension>::globalNumNodes() const {
  return allReduce(this->numNodes(), SPHERAL_OP_SUM);
}

//------------------------------------------------------------------------------
// Global numbers of fluid nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
size_t
DataBase<Dimension>::globalNumFluidInternalNodes() const {
  return allReduce(this->numFluidInternalNodes(), SPHERAL_OP_SUM);
}

template<typename Dimension>
size_t
DataBase<Dimension>::globalNumFluidGhostNodes() const {
  return allReduce(this->numFluidGhostNodes(), SPHERAL_OP_SUM);
}

template<typename Dimension>
size_t
DataBase<Dimension>::globalNumFluidNodes() const {
  return allReduce(this->numFluidNodes(), SPHERAL_OP_SUM);
}

//------------------------------------------------------------------------------
// Global numbers of solid nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
size_t
DataBase<Dimension>::globalNumSolidInternalNodes() const {
  return allReduce(this->numSolidInternalNodes(), SPHERAL_OP_SUM);
}

template<typename Dimension>
size_t
DataBase<Dimension>::globalNumSolidGhostNodes() const {
  return allReduce(this->numSolidGhostNodes(), SPHERAL_OP_SUM);
}

template<typename Dimension>
size_t
DataBase<Dimension>::globalNumSolidNodes() const {
  return allReduce(this->numSolidNodes(), SPHERAL_OP_SUM);
}

//------------------------------------------------------------------------------
// Global numbers of DEM nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
size_t
DataBase<Dimension>::globalNumDEMInternalNodes() const {
  return allReduce(this->numDEMInternalNodes(), SPHERAL_OP_SUM);
}

template<typename Dimension>
size_t
DataBase<Dimension>::globalNumDEMGhostNodes() const {
  return allReduce(this->numDEMGhostNodes(), SPHERAL_OP_SUM);
}

template<typename Dimension>
size_t
DataBase<Dimension>::globalNumDEMNodes() const {
  return allReduce(this->numDEMNodes(), SPHERAL_OP_SUM);
}

//------------------------------------------------------------------------------
// Optimize neighbor objects.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
reinitializeNeighbors() const {

  // // Find the current bounding box and average node extent in one loop.
  // // Compute the average node extent.
  // auto xmin = Vector(0.1*std::numeric_limits<Scalar>::max()), xmax = 0.1*Vector(std::numeric_limits<Scalar>::lowest());
  // unsigned ntot = 0;
  // Scalar havg = 0.0, hmax = 0.0, maxExtent = 0.0;
  // for (auto itr = this->nodeListBegin(); itr != this->nodeListEnd(); ++itr) {
  //   const auto positions = (*itr)->positions();
  //   const auto  n = (*itr)->numInternalNodes();
  //   const auto& neighbor = (*itr)->neighbor();
  //   const auto  etaMax = neighbor.kernelExtent();
  //   maxExtent = std::max(maxExtent, etaMax);
  //   ntot += n;
  //   for (auto i = 0u; i < n; ++i) {
  //     const auto& xi = positions(i);
  //     xmin = elementWiseMin(xmin, xi);
  //     xmax = elementWiseMax(xmax, xi);
  //     const auto hi = neighbor.nodeExtent(i).maxElement()/etaMax;
  //     havg += hi;
  //     hmax = std::max(hmax, hi);
  //   }
  // }

  // // Find the global result across all processors.
  // auto box = 0.0;
  // for (auto i = 0; i != Dimension::nDim; ++i) {
  //   xmin(i) = allReduce(xmin(i), SPHERAL_OP_MIN);
  //   xmax(i) = allReduce(xmax(i), SPHERAL_OP_MAX);
  //   box = std::max(box, xmax(i) - xmin(i));
  // }
  // havg = allReduce(havg, SPHERAL_OP_SUM);
  // ntot = allReduce(ntot, SPHERAL_OP_SUM);
  // hmax = allReduce(hmax, SPHERAL_OP_MAX);
  // if (ntot > 0) {
  //   havg /= ntot;

  //   box = std::max(box, maxExtent*hmax);
  //   for (auto i = 0; i < Dimension::nDim; ++i) {
  //     xmin(i) -= box;
  //     xmax(i) += box;
  //   }

  //   // Now initialize the neighbors.
  //   for (auto itr = this->nodeListBegin(); itr != this->nodeListEnd(); ++itr) {
  //     auto& neighbor = (*itr)->neighbor();
  //     neighbor.reinitialize(xmin, xmax, havg);
  //     neighbor.updateNodes();
  //   }
  // }
  VERIFY(false);
}

//------------------------------------------------------------------------------
// Update the connectivity map.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
updateConnectivityMap(const bool computeGhostConnectivity,
                      const bool computeOverlapConnectivity,
                      const bool computeIntersectionConnectivity) const {
  REQUIRE(mConnectivityMapPtr);
  mConnectivityMapPtr->rebuild(nodeListBegin(), nodeListEnd(),
                               computeGhostConnectivity, computeOverlapConnectivity, computeIntersectionConnectivity);
}

//------------------------------------------------------------------------------
// Patch the connectivity map.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
patchConnectivityMap(const FieldList<Dimension, size_t>& flags,
                     const FieldList<Dimension, size_t>& old2new) const {
  REQUIRE(mConnectivityMapPtr);
  mConnectivityMapPtr->patchConnectivity(flags, old2new);
}

//------------------------------------------------------------------------------
// Add a NodeList to this DataBase.
//------------------------------------------------------------------------------
// DEMNodeList
template<typename Dimension>
void
DataBase<Dimension>::
appendNodeList(DEMNodeList<Dimension>& nodeList) {
  REQUIRE(valid());
  if (this->haveNodeList(nodeList)) {
    SpheralWarning << "DataBase::appendNodeList() Warning: attempt to add DEMNodeList "
                   << &nodeList << " to DataBase " << this
                   << ", which already has it." << endl;
  } else {
    const NodeListRegistrar<Dimension>& nlr = NodeListRegistrar<Dimension>::instance();
    NodeListIterator orderItr = nlr.findInsertionPoint((NodeList<Dimension>*) &nodeList,
                                                       nodeListBegin(),
                                                       nodeListEnd());
    mNodeListPtrs.insert(orderItr, &nodeList);
    DEMNodeListIterator DEMOrderItr = nlr.findInsertionPoint(&nodeList,
                                                                 DEMNodeListBegin(),
                                                                 DEMNodeListEnd());
    const size_t delta = distance(DEMNodeListBegin(), DEMOrderItr);
    mDEMNodeListPtrs.insert(DEMOrderItr, &nodeList);
    mDEMNodeListAsNodeListPtrs.insert(DEMNodeListAsNodeListBegin() + delta,
                                        &nodeList);

  }
  updateView();
  ENSURE(valid());
}

// SolidNodeList
template<typename Dimension>
void
DataBase<Dimension>::
appendNodeList(SolidNodeList<Dimension>& nodeList) {
  REQUIRE(valid());
  if (this->haveNodeList(nodeList)) {
    SpheralWarning << "DataBase::appendNodeList() Warning: attempt to add SolidNodeList "
                   << &nodeList << " to DataBase " << this
                   << ", which already has it." << endl;
  } else {
    const NodeListRegistrar<Dimension>& nlr = NodeListRegistrar<Dimension>::instance();
    NodeListIterator orderItr = nlr.findInsertionPoint((NodeList<Dimension>*) &nodeList,
                                                       nodeListBegin(),
                                                       nodeListEnd());
    mNodeListPtrs.insert(orderItr, &nodeList);

    SolidNodeListIterator solidOrderItr = nlr.findInsertionPoint(&nodeList,
                                                                 solidNodeListBegin(),
                                                                 solidNodeListEnd());
    size_t delta = distance(solidNodeListBegin(), solidOrderItr);
    mSolidNodeListPtrs.insert(solidOrderItr, &nodeList);
    mSolidNodeListAsNodeListPtrs.insert(solidNodeListAsNodeListBegin() + delta,
                                        &nodeList);

    FluidNodeListIterator fluidOrderItr = nlr.findInsertionPoint(&nodeList,
                                                                 fluidNodeListBegin(),
                                                                 fluidNodeListEnd());
    delta = distance(fluidNodeListBegin(), fluidOrderItr);
    mFluidNodeListPtrs.insert(fluidOrderItr, &nodeList);
    mFluidNodeListAsNodeListPtrs.insert(fluidNodeListAsNodeListBegin() + delta,
                                        &nodeList);

  }
  updateView();
  ENSURE(valid());
}

// FluidNodeList
template<typename Dimension>
void
DataBase<Dimension>::
appendNodeList(FluidNodeList<Dimension>& nodeList) {
  REQUIRE(valid());
  if (this->haveNodeList(nodeList)) {
    SpheralWarning << "DataBase::appendNodeList() Warning: attempt to add FluidNodeList "
                   << &nodeList << " to DataBase " << this
                   << ", which already has it." << endl;
  } else {
    const NodeListRegistrar<Dimension>& nlr = NodeListRegistrar<Dimension>::instance();
    NodeListIterator orderItr = nlr.findInsertionPoint((NodeList<Dimension>*) &nodeList,
                                                       nodeListBegin(),
                                                       nodeListEnd());
    mNodeListPtrs.insert(orderItr, &nodeList);
    FluidNodeListIterator fluidOrderItr = nlr.findInsertionPoint(&nodeList,
                                                                 fluidNodeListBegin(),
                                                                 fluidNodeListEnd());
    const size_t delta = distance(fluidNodeListBegin(), fluidOrderItr);
    mFluidNodeListPtrs.insert(fluidOrderItr, &nodeList);
    mFluidNodeListAsNodeListPtrs.insert(fluidNodeListAsNodeListBegin() + delta,
                                        &nodeList);

  }
  updateView();
  ENSURE(valid());
}


// NodeList
template<typename Dimension>
void
DataBase<Dimension>::
appendNodeList(NodeList<Dimension>& nodeList) {
  REQUIRE(valid());
  if (this->haveNodeList(nodeList)) {
    SpheralWarning << "DataBase::appendNodeList() Warning: attempt to add NodeList "
                   << &nodeList << " to DataBase " << this
                   << ", which already has it." << endl;
  } else {
    const NodeListRegistrar<Dimension>& nlr = NodeListRegistrar<Dimension>::instance();
    NodeListIterator orderItr = nlr.findInsertionPoint(&nodeList,
                                                       nodeListBegin(),
                                                       nodeListEnd());
    mNodeListPtrs.insert(orderItr, &nodeList);
  }
  updateView();
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Delete a NodeList from this DataBase.
//------------------------------------------------------------------------------
// DEMNodeList
template<typename Dimension>
void
DataBase<Dimension>::
deleteNodeList(DEMNodeList<Dimension>& nodeList) {
  REQUIRE(valid());
  if (!this->haveNodeList(nodeList)) {
    SpheralWarning << "DataBase::deleteNodeList() Warning: attempt to remove DEMNodeList "
                   << &nodeList << " from DataBase " << this
                   << ", which does not have it." << endl;
  } else {
    // Erase from the NodeList vector.
    NodeListIterator nodeListItr = find(nodeListBegin(), nodeListEnd(),
                                        &nodeList);
    CHECK(nodeListItr != nodeListEnd());
    mNodeListPtrs.erase(nodeListItr);

    // Erase from the DEMNodeList vector.
    DEMNodeListIterator DEMItr = find(DEMNodeListBegin(), DEMNodeListEnd(), &nodeList);
    CHECK(DEMItr != DEMNodeListEnd());
    mDEMNodeListPtrs.erase(DEMItr);

    // Erase from the DEMNodeListAsNodeList vector.
    nodeListItr = find(DEMNodeListAsNodeListBegin(),
                       DEMNodeListAsNodeListEnd(),
                       &nodeList);
    CHECK(nodeListItr != DEMNodeListAsNodeListEnd());
    mDEMNodeListAsNodeListPtrs.erase(nodeListItr);
  }
  updateView();
  ENSURE(valid());
}

// SolidNodeList
template<typename Dimension>
void
DataBase<Dimension>::
deleteNodeList(SolidNodeList<Dimension>& nodeList) {
  REQUIRE(valid());
  if (!this->haveNodeList(nodeList)) {
    SpheralWarning << "DataBase::deleteNodeList() Warning: attempt to remove SolidNodeList "
                   << &nodeList << " from DataBase " << this
                   << ", which does not have it." << endl;
  } else {
    // Erase from the NodeList vector.
    NodeListIterator nodeListItr = find(nodeListBegin(), nodeListEnd(),
                                        &nodeList);
    CHECK(nodeListItr != nodeListEnd());
    mNodeListPtrs.erase(nodeListItr);

    // Erase from the SolidNodeList vector.
    SolidNodeListIterator solidItr = find(solidNodeListBegin(), solidNodeListEnd(), &nodeList);
    CHECK(solidItr != solidNodeListEnd());
    mSolidNodeListPtrs.erase(solidItr);

    // Erase from the SolidNodeListAsNodeList vector.
    nodeListItr = find(solidNodeListAsNodeListBegin(),
                       solidNodeListAsNodeListEnd(),
                       &nodeList);
    CHECK(nodeListItr != solidNodeListAsNodeListEnd());
    mSolidNodeListAsNodeListPtrs.erase(nodeListItr);

    // Erase from the FluidNodeList vector.
    FluidNodeListIterator fluidItr = find(fluidNodeListBegin(), fluidNodeListEnd(), &nodeList);
    CHECK(fluidItr != fluidNodeListEnd());
    mFluidNodeListPtrs.erase(fluidItr);

    // Erase from the FluidNodeListAsNodeList vector.
    nodeListItr = find(fluidNodeListAsNodeListBegin(),
                       fluidNodeListAsNodeListEnd(),
                       &nodeList);
    CHECK(nodeListItr != fluidNodeListAsNodeListEnd());
    mFluidNodeListAsNodeListPtrs.erase(nodeListItr);
  }
  updateView();
  ENSURE(valid());
}

// FluidNodeList
template<typename Dimension>
void
DataBase<Dimension>::
deleteNodeList(FluidNodeList<Dimension>& nodeList) {
  REQUIRE(valid());
  if (!this->haveNodeList(nodeList)) {
    SpheralWarning << "DataBase::deleteNodeList() Warning: attempt to remove FluidNodeList "
                   << &nodeList << " from DataBase " << this
                   << ", which does not have it." << endl;
  } else {
    // Erase from the NodeList vector.
    NodeListIterator nodeListItr = find(nodeListBegin(), nodeListEnd(),
                                        &nodeList);
    CHECK(nodeListItr != nodeListEnd());
    mNodeListPtrs.erase(nodeListItr);

    // Erase from the FluidNodeList vector.
    FluidNodeListIterator fluidItr = find(fluidNodeListBegin(), fluidNodeListEnd(), &nodeList);
    CHECK(fluidItr != fluidNodeListEnd());
    mFluidNodeListPtrs.erase(fluidItr);

    // Erase from the FluidNodeListAsNodeList vector.
    nodeListItr = find(fluidNodeListAsNodeListBegin(),
                       fluidNodeListAsNodeListEnd(),
                       &nodeList);
    CHECK(nodeListItr != fluidNodeListAsNodeListEnd());
    mFluidNodeListAsNodeListPtrs.erase(nodeListItr);
  }
  updateView();
  ENSURE(valid());
}

// NodeList
template<typename Dimension>
void
DataBase<Dimension>::
deleteNodeList(NodeList<Dimension>& nodeList) {
  REQUIRE(valid());
  if (!this->haveNodeList(nodeList)) {
    SpheralWarning << "DataBase::deleteNodeList() Warning: attempt to remove NodeList "
                   << &nodeList << " from DataBase " << this
                   << ", which does not have it." << endl;
  } else {
    // Erase from the NodeList vector.
    NodeListIterator nodeListItr = find(nodeListBegin(), nodeListEnd(),
                                        &nodeList);
    CHECK(nodeListItr != nodeListEnd());
    mNodeListPtrs.erase(nodeListItr);
  }
  updateView();
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Return the global mass field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::globalMass() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  for (auto* nptr: mNodeListPtrs) result.appendField(nptr->mass());
  return result;
}

//------------------------------------------------------------------------------
// Return the global mass (RZ) field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::globalMassRZ() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  if (GeometryRegistrar::coords() != CoordinateType::RZ) {
    SpheralWarning << "DataBase::globalMassRZ not valid in coordinate system " << GeometryRegistrar::coords() << std::endl;
  } else {
    for (auto* nptr: mNodeListPtrs) result.appendField(nptr->massRZ());
  }
  return result;
}

//------------------------------------------------------------------------------
// Return the global position field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
DataBase<Dimension>::globalPosition() const {
  REQUIRE(valid());
  FieldList<Dimension, Vector> result;
  for (auto* nptr: mNodeListPtrs) result.appendField(nptr->positions());
  return result;
}

//------------------------------------------------------------------------------
// Return the global velocity field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
DataBase<Dimension>::globalVelocity() const {
  REQUIRE(valid());
  FieldList<Dimension, Vector> result;
  for (auto* nptr: mNodeListPtrs) result.appendField(nptr->velocity());
  return result;
}

//------------------------------------------------------------------------------
// Return the global H smoothing field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::SymTensor>
DataBase<Dimension>::globalHfield() const {
  REQUIRE(valid());
  FieldList<Dimension, SymTensor> result;
  for (auto* nptr: mNodeListPtrs) result.appendField(nptr->Hfield());
  return result;
}

//------------------------------------------------------------------------------
// Return the global work field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::globalWork() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  for (auto* nptr: mNodeListPtrs) result.appendField(nptr->work());
  return result;
}

//------------------------------------------------------------------------------
// Return the fluid mass field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::fluidMass() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  for (auto* nptr: mFluidNodeListPtrs) result.appendField(nptr->mass());
  return result;
}

//------------------------------------------------------------------------------
// Return the fluid mass (RZ) field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::fluidMassRZ() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  if (GeometryRegistrar::coords() != CoordinateType::RZ) {
    SpheralWarning << "DataBase::fluidMassRZ not valid in coordinate system " << GeometryRegistrar::coords() << std::endl;
  } else {
    for (auto* nptr: mFluidNodeListPtrs) result.appendField(nptr->massRZ());
  }
  return result;
}

//------------------------------------------------------------------------------
// Return the fluid position field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
DataBase<Dimension>::fluidPosition() const {
  REQUIRE(valid());
  FieldList<Dimension, Vector> result;
  for (auto* nptr: mFluidNodeListPtrs) result.appendField(nptr->positions());
  return result;
}

//------------------------------------------------------------------------------
// Return the fluid velocity field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
DataBase<Dimension>::fluidVelocity() const {
  REQUIRE(valid());
  FieldList<Dimension, Vector> result;
  for (auto* nptr: mFluidNodeListPtrs) result.appendField(nptr->velocity());
  return result;
}

//------------------------------------------------------------------------------
// Return the fluid mass density field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::fluidMassDensity() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  for (auto* nptr: mFluidNodeListPtrs) result.appendField(nptr->massDensity());
  return result;
}

//------------------------------------------------------------------------------
// Return the fluid mass density (RZ) field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::fluidMassDensityRZ() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  if (GeometryRegistrar::coords() != CoordinateType::RZ) {
    SpheralWarning << "DataBase::fludMassDensityRZ not valid in coordinate system " << GeometryRegistrar::coords() << std::endl;
  } else {
    for (auto* nptr: mFluidNodeListPtrs) result.appendField(nptr->massDensityRZ());
  }
  return result;
}

//------------------------------------------------------------------------------
// Return the fluid specific thermal energy field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::fluidSpecificThermalEnergy() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  for (auto* nptr: mFluidNodeListPtrs) result.appendField(nptr->specificThermalEnergy());
  return result;
}

//------------------------------------------------------------------------------
// Return the fluid H smoothing field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::SymTensor>
DataBase<Dimension>::fluidHfield() const {
  REQUIRE(valid());
  FieldList<Dimension, SymTensor> result;
  for (auto* nptr: mFluidNodeListPtrs) result.appendField(nptr->Hfield());
  return result;
}

//------------------------------------------------------------------------------
// Return the fluid work field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::fluidWork() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  for (auto* nptr: mFluidNodeListPtrs) result.appendField(nptr->work());
  return result;
}

//------------------------------------------------------------------------------
// Return the solid mass field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::solidMass() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  for (auto* nptr: mSolidNodeListPtrs) result.appendField(nptr->mass());
  return result;
}

//------------------------------------------------------------------------------
// Return the solid mass (RZ) field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::solidMassRZ() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  if (GeometryRegistrar::coords() != CoordinateType::RZ) {
    SpheralWarning << "DataBase::solidMassRZ not valid in coordinate system " << GeometryRegistrar::coords() << std::endl;
  } else {
    for (auto* nptr: mSolidNodeListPtrs) result.appendField(nptr->massRZ());
  }
  return result;
}

//------------------------------------------------------------------------------
// Return the solid position field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
DataBase<Dimension>::solidPosition() const {
  REQUIRE(valid());
  FieldList<Dimension, Vector> result;
  for (auto* nptr: mSolidNodeListPtrs) result.appendField(nptr->positions());
  return result;
}

//------------------------------------------------------------------------------
// Return the solid velocity field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
DataBase<Dimension>::solidVelocity() const {
  REQUIRE(valid());
  FieldList<Dimension, Vector> result;
  for (auto* nptr: mSolidNodeListPtrs) result.appendField(nptr->velocity());
  return result;
}

//------------------------------------------------------------------------------
// Return the solid mass density field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::solidMassDensity() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  for (auto* nptr: mSolidNodeListPtrs) result.appendField(nptr->massDensity());
  return result;
}

//------------------------------------------------------------------------------
// Return the solid mass density (RZ) field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::solidMassDensityRZ() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  if (GeometryRegistrar::coords() != CoordinateType::RZ) {
    SpheralWarning << "DataBase::solidMassRZ not valid in coordinate system " << GeometryRegistrar::coords() << std::endl;
  } else {
    for (auto* nptr: mSolidNodeListPtrs) result.appendField(nptr->massDensityRZ());
  }
  return result;
}

//------------------------------------------------------------------------------
// Return the solid specific thermal energy field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::solidSpecificThermalEnergy() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  for (auto* nptr: mSolidNodeListPtrs) result.appendField(nptr->specificThermalEnergy());
  return result;
}

//------------------------------------------------------------------------------
// Return the solid H smoothing field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::SymTensor>
DataBase<Dimension>::solidHfield() const {
  REQUIRE(valid());
  FieldList<Dimension, SymTensor> result;
  for (auto* nptr: mSolidNodeListPtrs) result.appendField(nptr->Hfield());
  return result;
}

//------------------------------------------------------------------------------
// Return the solid work field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::solidWork() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  for (auto* nptr: mSolidNodeListPtrs) result.appendField(nptr->work());
  return result;
}

//------------------------------------------------------------------------------
// Return the solid deviatoric stress field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::SymTensor>
DataBase<Dimension>::solidDeviatoricStress() const {
  REQUIRE(valid());
  FieldList<Dimension, SymTensor> result;
  for (auto* nptr: mSolidNodeListPtrs) result.appendField(nptr->deviatoricStress());
  return result;
}

//------------------------------------------------------------------------------
// Return the solid plastic strain field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::solidPlasticStrain() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  for (auto* nptr: mSolidNodeListPtrs) result.appendField(nptr->plasticStrain());
  return result;
}

//------------------------------------------------------------------------------
// Return the solid plastic strain rate field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::solidPlasticStrainRate() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  for (auto* nptr: mSolidNodeListPtrs) result.appendField(nptr->plasticStrainRate());
  return result;
}

//------------------------------------------------------------------------------
// Return the solid damage field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::SymTensor>
DataBase<Dimension>::solidDamage() const {
  REQUIRE(valid());
  FieldList<Dimension, SymTensor> result;
  for (auto* nptr: mSolidNodeListPtrs) result.appendField(nptr->damage());
  return result;
}

//------------------------------------------------------------------------------
// Return the solid fragment ID field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, int>
DataBase<Dimension>::solidFragmentIDs() const {
  REQUIRE(valid());
  FieldList<Dimension, int> result;
  for (auto* nptr: mSolidNodeListPtrs) result.appendField(nptr->fragmentIDs());
  return result;
}

//------------------------------------------------------------------------------
// Return the solid particle type field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, int>
DataBase<Dimension>::solidParticleTypes() const {
  REQUIRE(valid());
  FieldList<Dimension, int> result;
  for (auto* nptr: mSolidNodeListPtrs) result.appendField(nptr->particleTypes());
  return result;
}

//------------------------------------------------------------------------------
// Return the DEM mass field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::DEMMass() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  for (auto* nptr: mDEMNodeListPtrs) result.appendField(nptr->mass());
  return result;
}

//------------------------------------------------------------------------------
// Return the DEM position field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
DataBase<Dimension>::DEMPosition() const {
  REQUIRE(valid());
  FieldList<Dimension, Vector> result;
  for (auto* nptr: mDEMNodeListPtrs) result.appendField(nptr->positions());
  return result;
}

//------------------------------------------------------------------------------
// Return the DEM velocity field.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
DataBase<Dimension>::DEMVelocity() const {
  REQUIRE(valid());
  FieldList<Dimension, Vector> result;
  for (auto* nptr: mDEMNodeListPtrs) result.appendField(nptr->velocity());
  return result;
}


//------------------------------------------------------------------------------
// Return the DEM HField.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::SymTensor>
DataBase<Dimension>::DEMHfield() const {
  REQUIRE(valid());
  FieldList<Dimension, SymTensor> result;
  for (auto* nptr: mDEMNodeListPtrs) result.appendField(nptr->Hfield());
  return result;
}


//------------------------------------------------------------------------------
// Return the DEM Particle Radii
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
DataBase<Dimension>::DEMParticleRadius() const {
  REQUIRE(valid());
  FieldList<Dimension, Scalar> result;
  for (auto* nptr: mDEMNodeListPtrs) result.appendField(nptr->particleRadius());
  return result;
}


//------------------------------------------------------------------------------
// Return the DEM unique particle index
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, size_t>
DataBase<Dimension>::DEMUniqueIndex() const {
  REQUIRE(valid());
  FieldList<Dimension, size_t> result;
  for (auto* nptr: mDEMNodeListPtrs) result.appendField(nptr->uniqueIndex());
  return result;
}


//------------------------------------------------------------------------------
// Return the DEM composite particle index
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, int>
DataBase<Dimension>::DEMCompositeParticleIndex() const {
  REQUIRE(valid());
  FieldList<Dimension, int> result;
  for (auto* nptr: mDEMNodeListPtrs) result.appendField(nptr->compositeParticleIndex());
  return result;
}

//------------------------------------------------------------------------------
// Set the DEM particle global indices
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::setDEMUniqueIndices() {
  REQUIRE(valid());
  auto uniqueIndices = this->DEMUniqueIndex();
  uniqueIndices += globalNodeIDs<Dimension>(this->nodeListBegin(),this->nodeListEnd());
}


//------------------------------------------------------------------------------
// Return the node extent for each NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
DataBase<Dimension>::globalNodeExtent() const {
  VERIFY(false);
  // REQUIRE(valid());
  // FieldList<Dimension, Vector> result;
  // for (ConstNodeListIterator nodeListItr = nodeListBegin();
  //      nodeListItr < nodeListEnd(); ++nodeListItr) {
  //   result.appendField((*nodeListItr)->neighbor().nodeExtentField());
  // }
  // return result;
}

//------------------------------------------------------------------------------
// Return the node extent for each FluidNodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
DataBase<Dimension>::fluidNodeExtent() const {
  VERIFY(false);
  // REQUIRE(valid());
  // FieldList<Dimension, Vector> result;
  // for (ConstFluidNodeListIterator nodeListItr = fluidNodeListBegin();
  //      nodeListItr < fluidNodeListEnd(); ++nodeListItr) {
  //   result.appendField((*nodeListItr)->neighbor().nodeExtentField());
  // }
  // return result;
}

//------------------------------------------------------------------------------
// Return the node extent for each SolidNodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
DataBase<Dimension>::solidNodeExtent() const {
  VERIFY(false);
  // REQUIRE(valid());
  // FieldList<Dimension, Vector> result;
  // for (ConstSolidNodeListIterator nodeListItr = solidNodeListBegin();
  //      nodeListItr < solidNodeListEnd(); ++nodeListItr) {
  //   result.appendField((*nodeListItr)->neighbor().nodeExtentField());
  // }
  // return result;
}

//------------------------------------------------------------------------------
// Return the node extent for each DEMNodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
DataBase<Dimension>::DEMNodeExtent() const {
  VERIFY(false);
  // REQUIRE(valid());
  // FieldList<Dimension, Vector> result;
  // for (ConstDEMNodeListIterator nodeListItr = DEMNodeListBegin();
  //      nodeListItr < DEMNodeListEnd(); ++nodeListItr) {
  //   result.appendField((*nodeListItr)->neighbor().nodeExtentField());
  // }
  // return result;
}

//------------------------------------------------------------------------------
// Return the global Hinverse field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
globalHinverse(FieldList<Dimension, typename Dimension::SymTensor>& result) const {
  REQUIRE(valid());
  this->resizeGlobalFieldList(result, SymTensor::zero());
  for (auto [k, nptr]: enumerate(mNodeListPtrs)) nptr->Hinverse(*result[k]);
}

//------------------------------------------------------------------------------
// Return the fluid Hinverse field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
fluidHinverse(FieldList<Dimension, typename Dimension::SymTensor>& result) const {
  REQUIRE(valid());
  this->resizeFluidFieldList(result, SymTensor::zero());
  for (auto [k, nptr]: enumerate(mFluidNodeListPtrs)) nptr->Hinverse(*result[k]);
}

//------------------------------------------------------------------------------
// Return the fluid pressure field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
fluidPressure(FieldList<Dimension, typename Dimension::Scalar>& result) const {
  REQUIRE(valid());
  this->resizeFluidFieldList(result, 0.0, HydroFieldNames::pressure);
  for (auto [k, nptr]: enumerate(mFluidNodeListPtrs)) nptr->pressure(*result[k]);
}

//------------------------------------------------------------------------------
// Return the fluid temperature field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
fluidTemperature(FieldList<Dimension, typename Dimension::Scalar>& result) const {
  REQUIRE(valid());
  this->resizeFluidFieldList(result, 0.0, HydroFieldNames::temperature);
  for (auto [k, nptr]: enumerate(mFluidNodeListPtrs)) nptr->temperature(*result[k]);
}

//------------------------------------------------------------------------------
// Return the fluid sound speed field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
fluidSoundSpeed(FieldList<Dimension, typename Dimension::Scalar>& result) const {
  REQUIRE(valid());
  this->resizeFluidFieldList(result, 0.0, HydroFieldNames::soundSpeed);
  for (auto [k, nptr]: enumerate(mFluidNodeListPtrs)) nptr->soundSpeed(*result[k]);
}

//------------------------------------------------------------------------------
// Return the fluid volume field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
fluidVolume(FieldList<Dimension, typename Dimension::Scalar>& result) const {
  REQUIRE(valid());
  this->resizeFluidFieldList(result, 0.0, HydroFieldNames::volume, false);
  for (auto [k, nptr]: enumerate(mFluidNodeListPtrs)) nptr->volume(*result[k]);
}

//------------------------------------------------------------------------------
// Return the fluid gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
fluidGamma(FieldList<Dimension, typename Dimension::Scalar>& result) const {
  REQUIRE(valid());
  this->resizeFluidFieldList(result, 0.0, HydroFieldNames::gamma, false);
  for (auto [k, nptr]: enumerate(mFluidNodeListPtrs)) {
    const auto& eos = nptr->equationOfState();
    const auto& rho = nptr->massDensity();
    const auto& eps = nptr->specificThermalEnergy();
    eos.setGammaField(*result[k], rho, eps);
  }
}

//------------------------------------------------------------------------------
// Return the fluid entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
fluidEntropy(FieldList<Dimension, typename Dimension::Scalar>& result) const {
  REQUIRE(valid());
  this->resizeFluidFieldList(result, 0.0, HydroFieldNames::entropy, false);
  for (auto [k, nptr]: enumerate(mFluidNodeListPtrs)) {
    const auto& eos = nptr->equationOfState();
    const auto& rho = nptr->massDensity();
    const auto& eps = nptr->specificThermalEnergy();
    eos.setEntropy(*result[k], rho, eps);
  }
}

//------------------------------------------------------------------------------
// Return the fluid linear momentum field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
fluidLinearMomentum(FieldList<Dimension, typename Dimension::Vector>& result) const {
  REQUIRE(valid());
  this->resizeFluidFieldList(result, Vector::zero(), HydroFieldNames::linearMomentum);
  for (auto [k, nptr]: enumerate(mFluidNodeListPtrs)) nptr->linearMomentum(*result[k]);
}

//------------------------------------------------------------------------------
// Return the fluid energy field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
fluidTotalEnergy(FieldList<Dimension, typename Dimension::Scalar>& result) const {
  REQUIRE(valid());
  this->resizeFluidFieldList(result, 0.0, HydroFieldNames::totalEnergy, false);
  for (auto [k, nptr]: enumerate(mFluidNodeListPtrs)) nptr->totalEnergy(*result[k]);
}

//------------------------------------------------------------------------------
// Return the specific heat field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
fluidSpecificHeat(const FieldList<Dimension, typename Dimension::Scalar>& temperature,
                  FieldList<Dimension, typename Dimension::Scalar>& result) const {
  REQUIRE(valid());
  this->resizeFluidFieldList(result, 0.0, HydroFieldNames::specificHeat);
  for (auto [k, nptr]: enumerate(mFluidNodeListPtrs)) {
    const auto& eos = nptr->equationOfState();
    const auto& rho = nptr->massDensity();
    const auto& temp = *temperature[k];
    eos.setSpecificHeat(*result[k], rho, temp);
  }
}

//------------------------------------------------------------------------------
// Return the number of neighbors for each node based on the last update of 
// the connectivity map.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, size_t>
DataBase<Dimension>::numNeighbors() const {
  REQUIRE(valid());
  VERIFY(mConnectivityMapPtr);
  FieldList<Dimension, size_t> result = newFluidFieldList(size_t(0u), "number of neighbors");
  const auto& pairs = mConnectivityMapPtr->nodePairList();
  for (const auto& p: pairs) {
    ++result(p.i_list, p.i_node);
    ++result(p.j_list, p.j_node);
  }
  return result;
}

//------------------------------------------------------------------------------
// Return the number of internal nodes in sets of NodeLists
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<size_t>
DataBase<Dimension>::numInternalNodesPerNodeList() const {
  const auto n = mNodeListPtrs.size();
  std::vector<size_t> result(n);
  for (auto i = 0u; i < n; ++i) result[i] = mNodeListPtrs[i]->numInternalNodes();
  return result;
}

template<typename Dimension>
std::vector<size_t>
DataBase<Dimension>::numInternalNodesPerFluidNodeList() const {
  const auto n = mFluidNodeListPtrs.size();
  std::vector<size_t> result(n);
  for (auto i = 0u; i < n; ++i) result[i] = mFluidNodeListPtrs[i]->numInternalNodes();
  return result;
}

template<typename Dimension>
std::vector<size_t>
DataBase<Dimension>::numInternalNodesPerSolidNodeList() const {
  const auto n = mSolidNodeListPtrs.size();
  std::vector<size_t> result(n);
  for (auto i = 0u; i < n; ++i) result[i] = mSolidNodeListPtrs[i]->numInternalNodes();
  return result;
}

template<typename Dimension>
std::vector<size_t>
DataBase<Dimension>::numInternalNodesPerDEMNodeList() const {
  const auto n = mDEMNodeListPtrs.size();
  std::vector<size_t> result(n);
  for (auto i = 0u; i < n; ++i) result[i] = mDEMNodeListPtrs[i]->numInternalNodes();
  return result;
}

//------------------------------------------------------------------------------
// Compute coordinates bounding all nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::boundingBox(Vector& xmin,
                                 Vector& xmax,
                                 const bool ghost) const {

  // This method is now exported to an external stand-alone.  Maintained here
  // for backwards compatibility.
  const FieldList<Dimension, Vector> positions = this->globalPosition();
  Spheral::globalBoundingBox(positions, xmin, xmax, ghost);
}

//------------------------------------------------------------------------------
// Compute coordinates bounding all nodes in the DataBase.
// This version allows the user to mask off nodes (mask = 0 to ignore).
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::boundingBox(Vector& xmin,
                                 Vector& xmax,
                                 const FieldList<Dimension, int>& mask,
                                 const bool ghost) const {

  // Find our local bounds.
  xmin = DBL_MAX;
  xmax = -DBL_MAX;
  const FieldList<Dimension, Vector> positions = this->globalPosition();
  for (auto nodeList = 0u; nodeList < positions.numFields(); ++nodeList) {
    const auto n = ghost ? positions[nodeList]->numElements() : positions[nodeList]->numInternalElements();
    for (auto i = 0u; i < n; ++i) {
      if (mask(nodeList, i) != 0) {
        const auto& xi = positions(nodeList, i);
        xmin = elementWiseMin(xmin, xi);
        xmax = elementWiseMax(xmax, xi);
      }
    }
  }

  // Now find the global bounds across all processors.
  for (auto i = 0; i < Dimension::nDim; ++i) {
    xmin(i) = allReduce(xmin(i), SPHERAL_OP_MIN);
    xmax(i) = allReduce(xmax(i), SPHERAL_OP_MAX);
  }
}

//------------------------------------------------------------------------------
// Return the local (per domain) sampling extents.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::localSamplingBoundingVolume(Vector& centroid,
                                                 double& radiusNodes,
                                                 double& radiusSample,
                                                 Vector& xminNodes,
                                                 Vector& xmaxNodes,
                                                 Vector& xminSample,
                                                 Vector& xmaxSample) const {

  // Find our local bounds.
  centroid.Zero();
  xminNodes = FLT_MAX;
  xmaxNodes = -FLT_MAX;
  xminSample = FLT_MAX;
  xmaxSample = -FLT_MAX;
  size_t count = 0;
  const FieldList<Dimension, Vector> positions = this->globalPosition();
  const FieldList<Dimension, Vector> extent = this->globalNodeExtent();
  for (auto nodeList = 0u; nodeList != positions.numFields(); ++nodeList) {
    for (auto i = 0u; i != mNodeListPtrs[nodeList]->numInternalNodes(); ++i) {
      const Vector& xi = positions(nodeList, i);
      const Vector& extenti = extent(nodeList, i);
      const Vector xmini = xi - extenti;
      const Vector xmaxi = xi + extenti;
      centroid += xi;
      xminNodes = elementWiseMin(xminNodes, xi);
      xmaxNodes = elementWiseMax(xmaxNodes, xi);
      xminSample = elementWiseMin(xminSample, xmini);
      xmaxSample = elementWiseMax(xmaxSample, xmaxi);
      ++count;
    }
  }

  // Normalize the centroid.
  if (count > 0) centroid /= count;

  // Find the maximal radial extent from the centroid.
  radiusNodes = 0.0;
  radiusSample = 0.0;
  FieldList<Dimension, SymTensor> Hinv = newGlobalFieldList(SymTensor::zero());
  this->globalHinverse(Hinv);
  for (auto nodeList = 0u; nodeList != positions.numFields(); ++nodeList) {
    for (auto i = 0u; i != mNodeListPtrs[nodeList]->numInternalNodes(); ++i) {
      const Vector& xi = positions(nodeList, i);
      const Vector dr = xi - centroid;
      const double drMag = dr.magnitude();
      const double hi = (Hinv(nodeList, i)*dr).magnitude() * safeInv(drMag, 1.0e-20);
      radiusNodes = std::max(radiusNodes, drMag);
      radiusSample = std::max(radiusSample, drMag + 2.0*hi);
    }
  }

  // Puff things out a tiny bit in an effort to make sure we're domain 
  // decomposition independent.
  const Vector delta = 0.001*(xmaxSample - xminSample);
  xminNodes -= delta;
  xmaxNodes += delta;
  xminSample -= delta;
  xmaxSample += delta;
  radiusNodes *= 1.001;
  radiusSample *= 1.001;
}

//------------------------------------------------------------------------------
// Return the global min and max sampling extents.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
globalSamplingBoundingVolume(Vector& centroid,
                             double& radiusNodes,
                             double& radiusSample,
                             Vector& xminNodes,
                             Vector& xmaxNodes,
                             Vector& xminSample,
                             Vector& xmaxSample) const {

  // Find our local bounds.
  this->localSamplingBoundingVolume(centroid, radiusNodes, radiusSample, 
                                    xminNodes, xmaxNodes,
                                    xminSample, xmaxSample);

#ifdef SPHERAL_ENABLE_MPI
  // Now find the global bounds across all processors.
  {
    size_t nlocal = this->numInternalNodes();
    centroid *= nlocal;
    for (int i = 0; i != Dimension::nDim; ++i) {
      xminNodes(i) = allReduce(xminNodes(i), SPHERAL_OP_MIN);
      xmaxNodes(i) = allReduce(xmaxNodes(i), SPHERAL_OP_MAX);
      xminSample(i) = allReduce(xminSample(i), SPHERAL_OP_MIN);
      xmaxSample(i) = allReduce(xmaxSample(i), SPHERAL_OP_MAX);
      centroid(i) = allReduce(centroid(i), SPHERAL_OP_SUM);
    }

    // Fix up the centroid and radii.
    size_t nglobal = allReduce((uint64_t) nlocal, SPHERAL_OP_SUM);
    if (nglobal > 0) {
      centroid /= nglobal;
      radiusNodes = 0.0;
      radiusSample = 0.0;
      const FieldList<Dimension, Vector> positions = this->globalPosition();
      const FieldList<Dimension, Vector> extent = this->globalNodeExtent();
      FieldList<Dimension, SymTensor> Hinv = this->newGlobalFieldList(SymTensor::zero(), "H inverse");
      this->globalHinverse(Hinv);
      for (auto nodeList = 0u; nodeList != positions.numFields(); ++nodeList) {
        for (auto i = 0u; i != mNodeListPtrs[nodeList]->numInternalNodes(); ++i) {
          const Vector& xi = positions(nodeList, i);
          const Vector dr = xi - centroid;
          const Vector drUnit = dr.unitVector();
          const double drMag = dr.magnitude();
          const double hi = (Hinv(nodeList, i)*drUnit).magnitude();
          radiusNodes = std::max(radiusNodes, drMag);
          radiusSample = std::max(radiusSample, drMag + 2.0*hi);
        }
      }
      radiusNodes = allReduce(radiusNodes, SPHERAL_OP_MAX);
      radiusSample = allReduce(radiusSample, SPHERAL_OP_MAX);
      radiusNodes *= 1.001;
      radiusSample *= 1.001;
    }
  }
#endif
}

//------------------------------------------------------------------------------
// Return the local (per domain) min and max sampling extents for nodes that
// are connected.  This could represent an arbitrary number of results (in the
// extreme case where all nodes a disconnected you get back a set the size of 
// the number of nodes!)
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::localSamplingBoundingBoxes(vector<Vector>& xminima,
                                                vector<Vector>& xmaxima) const {

  // Initialize the result.
  xminima = vector<Vector>();
  xmaxima = vector<Vector>();

  // We use our connectivity to make this more efficient.
  this->updateConnectivityMap(false, false, false);
  const ConnectivityMap<Dimension>& connectivityMap = this->connectivityMap(false, false, false);
  const FieldList<Dimension, Vector> positions = this->globalPosition();
  const FieldList<Dimension, Vector> extent = this->globalNodeExtent();
  const int numNodeLists = this->numNodeLists();

  // Iterate over the nodes and progressively build up the result by linking
  // things together as they're connected.
  FieldList<Dimension, int> flags = this->newGlobalFieldList(0, "flags");
  for (int nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = positions[nodeListi]->nodeList();
    for (auto i = 0u; i != nodeList.numInternalNodes(); ++i) {
      if (flags(nodeListi, i) == 0) {
        const Vector& xi = positions(nodeListi, i);
        const Vector& extenti = extent(nodeListi, i);
        Vector xminGroup = xi - extenti;
        Vector xmaxGroup = xi + extenti;
        flags(nodeListi, i) = 1;

        // Find the min/max extent of this connected set of nodes.
        const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);
        CHECK((int)fullConnectivity.size() == numNodeLists);
        for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
          const vector<int>& connectivity = fullConnectivity[nodeListj];
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;
            const Vector& xj = positions(nodeListj, j);
            const Vector& extentj = extent(nodeListj, j);
            xminGroup = elementWiseMin(xminGroup, xj - extentj);
            xmaxGroup = elementWiseMax(xmaxGroup, xj + extentj);
            flags(nodeListj, j) = 1;
          }
        }

        // Compare this batches min/max extents with the results so far.
        CHECK(xminima.size() == xmaxima.size());
        unsigned int k = 0;
        bool matched = false;
        while (k != xminima.size() and not matched) {
          if (testBoxIntersection(xminGroup, xmaxGroup, xminima[k], xmaxima[k])) {
            xminima[k] = elementWiseMin(xminima[k], xminGroup);
            xmaxima[k] = elementWiseMax(xmaxima[k], xmaxGroup);
            matched = true;
          }
          ++k;
        }
        if (not matched) {
          xminima.push_back(xminGroup);
          xmaxima.push_back(xmaxGroup);
        }
      }
    }
  }

  // Look for any overlaps between the boxes.
  int ibox = xminima.size() - 1;
  while (ibox > 0) {
    int jbox = 0;
    while (jbox != ibox) {
      if (testBoxIntersection(xminima[ibox], xmaxima[ibox], xminima[jbox], xmaxima[jbox])) {
        xminima[jbox] = elementWiseMin(xminima[jbox], xminima[ibox]);
        xmaxima[jbox] = elementWiseMax(xmaxima[jbox], xmaxima[ibox]);
        xminima.erase(xminima.begin() + ibox);
        xmaxima.erase(xmaxima.begin() + ibox);
        jbox = ibox;
      } else {
        ++jbox;
      }
    }
    --ibox;
  }

}

//------------------------------------------------------------------------------
// The global version of above.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::
globalSamplingBoundingBoxes(vector<Vector>& xminima,
                            vector<Vector>& xmaxima) const {

  // First get each domains local values.
  localSamplingBoundingBoxes(xminima, xmaxima);

#ifdef SPHERAL_ENABLE_MPI
  // Parallel crap.
  const int procID = Process::getRank();
  const int numProcs = Process::getTotalNumberOfProcesses();

  // Get each domains min & maxes.
  vector<char> localBuffer;
  packElement(xminima, localBuffer);
  packElement(xmaxima, localBuffer);
  xminima = vector<Vector>();
  xmaxima = vector<Vector>();
  for (int sendProc = 0; sendProc != numProcs; ++sendProc) {
    vector<char> buffer = localBuffer;
    int bufSize = localBuffer.size();
    MPI_Bcast(&bufSize, 1, MPI_INT, sendProc, Communicator::communicator());
    if (procID != sendProc) buffer.resize(bufSize);
    MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, Communicator::communicator());
    vector<Vector> xminimaDomain, xmaximaDomain;
    vector<char>::const_iterator bufItr = buffer.begin();
    unpackElement(xminimaDomain, bufItr, buffer.end());
    unpackElement(xmaximaDomain, bufItr, buffer.end());
    CHECK(bufItr == buffer.end());
    copy(xminimaDomain.begin(), xminimaDomain.end(), back_inserter(xminima));
    copy(xmaximaDomain.begin(), xmaximaDomain.end(), back_inserter(xmaxima));
  }

  // Go through the boxes for each domain, joining them up where possible.
  int ibox = xminima.size() - 1;
  while (ibox > 0) {
    int jbox = 0;
    while (jbox != ibox) {
      if (testBoxIntersection(xminima[ibox], xmaxima[ibox], xminima[jbox], xmaxima[jbox])) {
        xminima[jbox] = elementWiseMin(xminima[jbox], xminima[ibox]);
        xmaxima[jbox] = elementWiseMax(xmaxima[jbox], xmaxima[ibox]);
        xminima.erase(xminima.begin() + ibox);
        xmaxima.erase(xmaxima.begin() + ibox);
        jbox = ibox;
      } else {
        ++jbox;
      }
    }
    --ibox;
  }
#endif
}

//------------------------------------------------------------------------------
// Test if the DataBase is "valid", or internally consistent.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
DataBase<Dimension>::valid() const {

  bool ok = (this->numNodeLists() >= this->numFluidNodeLists());

  // Verify that all FluidNodeLists are listed in the NodeList array.
  ConstFluidNodeListIterator fluidNodeItr = fluidNodeListBegin();
  while (ok && fluidNodeItr < fluidNodeListEnd()) {
    ok = this->haveNodeList(dynamic_cast<const NodeList<Dimension>&>(**fluidNodeItr));
    ++fluidNodeItr;
  }

  // Verify that all FluidNodeLists are listed in the FluidNodeListAsNodeList
  // array.
  fluidNodeItr = fluidNodeListBegin();
  while (ok && fluidNodeItr < fluidNodeListEnd()) {
    ok = find(mFluidNodeListAsNodeListPtrs.begin(),
              mFluidNodeListAsNodeListPtrs.end(),
              *fluidNodeItr) != mFluidNodeListAsNodeListPtrs.end();
    ++fluidNodeItr;
  }

  // Verify that all SolidNodeLists are listed in the NodeList array.
  ConstSolidNodeListIterator solidNodeItr = solidNodeListBegin();
  while (ok && solidNodeItr < solidNodeListEnd()) {
    ok = this->haveNodeList(dynamic_cast<const NodeList<Dimension>&>(**solidNodeItr));
    ++solidNodeItr;
  }

  // Verify that all SolidNodeLists are listed in the SolidNodeListAsNodeList
  // array.
  solidNodeItr = solidNodeListBegin();
  while (ok && solidNodeItr < solidNodeListEnd()) {
    ok = find(mSolidNodeListAsNodeListPtrs.begin(),
              mSolidNodeListAsNodeListPtrs.end(),
              *solidNodeItr) != mSolidNodeListAsNodeListPtrs.end();
    ++solidNodeItr;
  }

  // Verify that the NodeLists are listed in the proper order.
  const NodeListRegistrar<Dimension>& nlr = NodeListRegistrar<Dimension>::instance();
  if (ok) {
    typename NodeListRegistrar<Dimension>::const_iterator checkItr = nlr.begin();
    ConstNodeListIterator itr = nodeListBegin();
    while (itr != nodeListEnd()) {
      while (checkItr != nlr.end() && *checkItr != *itr)
        ++checkItr;
      ++itr;
    }
    ok = (itr == nodeListEnd());
  }

  // Verify that the FluidNodeLists are listed in the proper order.
  if (ok) {
    typename NodeListRegistrar<Dimension>::const_iterator checkItr = nlr.begin();
    ConstFluidNodeListIterator itr = fluidNodeListBegin();
    while (itr != fluidNodeListEnd()) {
      while (checkItr != nlr.end() && *checkItr != dynamic_cast<NodeList<Dimension>*>(*itr))
        ++checkItr;
      ++itr;
    }
    ok = (itr == fluidNodeListEnd());
  }

  // Verify that the FluidNodeListsAsNodeLists are listed in the proper order.
  if (ok) {
    typename NodeListRegistrar<Dimension>::const_iterator checkItr = nlr.begin();
    ConstNodeListIterator itr = fluidNodeListAsNodeListBegin();
    while (itr != fluidNodeListAsNodeListEnd()) {
      while (checkItr != nlr.end() && *checkItr != *itr)
        ++checkItr;
      ++itr;
    }
    ok = (itr == fluidNodeListAsNodeListEnd());
  }

  // Verify that the SolidNodeLists are listed in the proper order.
  if (ok) {
    typename NodeListRegistrar<Dimension>::const_iterator checkItr = nlr.begin();
    ConstSolidNodeListIterator itr = solidNodeListBegin();
    while (itr != solidNodeListEnd()) {
      while (checkItr != nlr.end() && *checkItr != dynamic_cast<NodeList<Dimension>*>(*itr))
        ++checkItr;
      ++itr;
    }
    ok = (itr == solidNodeListEnd());
  }

  // Verify that the SolidNodeListsAsNodeLists are listed in the proper order.
  if (ok) {
    typename NodeListRegistrar<Dimension>::const_iterator checkItr = nlr.begin();
    ConstNodeListIterator itr = solidNodeListAsNodeListBegin();
    while (itr != solidNodeListAsNodeListEnd()) {
      while (checkItr != nlr.end() && *checkItr != *itr)
        ++checkItr;
      ++itr;
    }
    ok = (itr == solidNodeListAsNodeListEnd());
  }

  // Check that the View containers are correct
  if (ok) {
    ok &= mNodeListViewStorage.size() == mNodeListPtrs.size();
    ok &= mFluidNodeListViewStorage.size() == mFluidNodeListPtrs.size();
    ok &= mSolidNodeListViewStorage.size() == mSolidNodeListPtrs.size();
    ok &= mDEMNodeListViewStorage.size() == mDEMNodeListPtrs.size();
    if (ok) {
      for (auto k = 0u; k < mNodeListViewStorage.size(); ++k)      ok &= mNodeListViewStorage[k] == static_cast<NodeListView<Dimension>&>(*mNodeListPtrs[k]);
      for (auto k = 0u; k < mFluidNodeListViewStorage.size(); ++k) ok &= mFluidNodeListViewStorage[k] == static_cast<FluidNodeListView<Dimension>&>(*mFluidNodeListPtrs[k]);
      for (auto k = 0u; k < mSolidNodeListViewStorage.size(); ++k) ok &= mSolidNodeListViewStorage[k] == static_cast<SolidNodeListView<Dimension>&>(*mSolidNodeListPtrs[k]);
      for (auto k = 0u; k < mDEMNodeListViewStorage.size(); ++k)   ok &= mDEMNodeListViewStorage[k] == static_cast<DEMNodeListView<Dimension>&>(*mDEMNodeListPtrs[k]);
    }
  }

  return ok;
}

//------------------------------------------------------------------------------
// Set the view containers
//------------------------------------------------------------------------------
template<typename Dimension>
void
DataBase<Dimension>::updateView() {

  // Clear out view storage
#ifdef SPHERAL_UNIFIED_MEMORY
  mNodeListViewStorage.resize(mNodeListPtrs.size());
  mFluidNodeListViewStorage.resize(mFluidNodeListPtrs.size());
  mSolidNodeListViewStorage.resize()mSolidNodeListPtrs.size();
  mDEMNodeListViewStorage.resize(mDEMNodeListPtrs.size());
#else
  auto allocateManagedArray = []<typename T>(ViewSpan<T>& array, const size_t n) {
    if (array.size() == 0u and !array.data(chai::CPU, false)) {
      array.allocate(n, chai::CPU);
    } else {
      array.reallocate(n);
    }
  };
  allocateManagedArray(mNodeListViewStorage,      mNodeListPtrs.size());
  allocateManagedArray(mFluidNodeListViewStorage, mFluidNodeListPtrs.size());
  allocateManagedArray(mSolidNodeListViewStorage, mSolidNodeListPtrs.size());
  allocateManagedArray(mDEMNodeListViewStorage,   mDEMNodeListPtrs.size());
#endif

  // Set storage arrays for View
  for (auto [i, nptr]: enumerate(mNodeListPtrs))      mNodeListViewStorage[i] = nptr->view();
  for (auto [i, nptr]: enumerate(mFluidNodeListPtrs)) mFluidNodeListViewStorage[i] = nptr->view();
  for (auto [i, nptr]: enumerate(mSolidNodeListPtrs)) mSolidNodeListViewStorage[i] = nptr->view();
  for (auto [i, nptr]: enumerate(mDEMNodeListPtrs))   mDEMNodeListViewStorage[i] = nptr->view();

  // Assign View spans
#ifdef SPHERAL_UNIFIED_MEMORY
  DataBaseView<Dimension>::mNodeListViews      = mNodeListViewStorage;
  DataBaseView<Dimension>::mFluidNodeListViews = mFluidNodeListViewStorage;
  DataBaseView<Dimension>::mSolidNodeListViews = mSolidNodeListViewStorage;
  DataBaseView<Dimension>::mDEMNodeListViews   = mDEMNodeListViewStorage;
#else
  for (auto i = 0u; i < mNodeListViewStorage.size(); ++i)      DataBaseView<Dimension>::mNodeListViews[i]      = mNodeListViewStorage[i];
  for (auto i = 0u; i < mFluidNodeListViewStorage.size(); ++i) DataBaseView<Dimension>::mFluidNodeListViews[i] = mFluidNodeListViewStorage[i];
  for (auto i = 0u; i < mSolidNodeListViewStorage.size(); ++i) DataBaseView<Dimension>::mSolidNodeListViews[i] = mSolidNodeListViewStorage[i];
  for (auto i = 0u; i < mDEMNodeListViewStorage.size(); ++i)   DataBaseView<Dimension>::mDEMNodeListViews[i]   = mDEMNodeListViewStorage[i];
#endif

  ENSURE(valid());
}

}

