//---------------------------------Spheral++----------------------------------//
// DataBaseView -- A device capable view of DataBase
//
// Created by JMO, Tue Oct 20 13:45:31 PDT 2026
//----------------------------------------------------------------------------//
#include "DataBase/DataBaseView.hh"
#include "Field/FieldView.hh"
#include "Field/FieldListView.hh"
#include "Utilities/DBC.hh"

#include "NodeList/NodeListView.hh"
#include "NodeList/FluidNodeListView.hh"
#include "NodeList/SolidNodeListView.hh"
#include "NodeList/DEMNodeListView.hh"

#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Numbers of nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::numInternalNodes() const {
  return std::accumulate(mNodeListViews.begin(), mNodeListViews.end(), size_t(0u),
                         [](size_t sum, const auto& x) { return sum + x.numInternalNodes(); });
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::numGhostNodes() const {
  return std::accumulate(mNodeListViews.begin(), mNodeListViews.end(), size_t(0u),
                         [](size_t sum, const auto& x) { return sum + x.numGhostNodes(); });
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::numNodes() const {
  return std::accumulate(mNodeListViews.begin(), mNodeListViews.end(), size_t(0u),
                         [](size_t sum, const auto& x) { return sum + x.numNodes(); });
}

//------------------------------------------------------------------------------
// Numbers of fluid nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::numFluidInternalNodes() const {
  return std::accumulate(mFluidNodeListViews.begin(), mFluidNodeListViews.end(), size_t(0u),
                         [](size_t sum, const auto& x) { return sum + x.numInternalNodes(); });
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::numFluidGhostNodes() const {
  return std::accumulate(mFluidNodeListViews.begin(), mFluidNodeListViews.end(), size_t(0u),
                         [](size_t sum, const auto& x) { return sum + x.numGhostNodes(); });
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::numFluidNodes() const {
  return std::accumulate(mFluidNodeListViews.begin(), mFluidNodeListViews.end(), size_t(0u),
                         [](size_t sum, const auto& x) { return sum + x.numNodes(); });
}

//------------------------------------------------------------------------------
// Numbers of solid nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::numSolidInternalNodes() const {
  return std::accumulate(mSolidNodeListViews.begin(), mSolidNodeListViews.end(), size_t(0u),
                         [](size_t sum, const auto& x) { return sum + x.numInternalNodes(); });
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::numSolidGhostNodes() const {
  return std::accumulate(mSolidNodeListViews.begin(), mSolidNodeListViews.end(), size_t(0u),
                         [](size_t sum, const auto& x) { return sum + x.numGhostNodes(); });
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::numSolidNodes() const {
  return std::accumulate(mSolidNodeListViews.begin(), mSolidNodeListViews.end(), size_t(0u),
                         [](size_t sum, const auto& x) { return sum + x.numNodes(); });
}

//------------------------------------------------------------------------------
// Numbers of DEM nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::numDEMInternalNodes() const {
  return std::accumulate(mDEMNodeListViews.begin(), mDEMNodeListViews.end(), size_t(0u),
                         [](size_t sum, const auto& x) { return sum + x.numInternalNodes(); });
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::numDEMGhostNodes() const {
  return std::accumulate(mDEMNodeListViews.begin(), mDEMNodeListViews.end(), size_t(0u),
                         [](size_t sum, const auto& x) { return sum + x.numGhostNodes(); });
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::numDEMNodes() const {
  return std::accumulate(mDEMNodeListViews.begin(), mDEMNodeListViews.end(), size_t(0u),
                         [](size_t sum, const auto& x) { return sum + x.numNodes(); });
}

//------------------------------------------------------------------------------
// Test if the given NodeList is registered with this DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
bool
DataBaseView<Dimension>::
haveNodeList(const NodeListView<Dimension>& x) const {
  return std::ranges::find(mNodeListViews, x) != mNodeListViews.end();
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
bool
DataBaseView<Dimension>::
haveNodeList(const FluidNodeListView<Dimension>& x) const {
  return std::ranges::find(mFluidNodeListViews, x) != mFluidNodeListViews.end();
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
bool
DataBaseView<Dimension>::
haveNodeList(const SolidNodeListView<Dimension>& x) const {
  return std::ranges::find(mSolidNodeListViews, x) != mSolidNodeListViews.end();
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
bool
DataBaseView<Dimension>::
haveNodeList(const DEMNodeListView<Dimension>& x) const {
  return std::ranges::find(mDEMNodeListViews, x) != mDEMNodeListViews.end();
}

//------------------------------------------------------------------------------
// Get the index for the given NodeList in the DataBase arrays
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::
nodeListIndex(const NodeListView<Dimension>& x) const {
  return std::distance(mNodeListViews.begin(), std::ranges::find(mNodeListViews, x));
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::
nodeListIndex(const FluidNodeListView<Dimension>& x) const {
  return std::distance(mFluidNodeListViews.begin(), std::ranges::find(mFluidNodeListViews, x));
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::
nodeListIndex(const SolidNodeListView<Dimension>& x) const {
  return std::distance(mSolidNodeListViews.begin(), std::ranges::find(mSolidNodeListViews, x));
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
DataBaseView<Dimension>::
nodeListIndex(const DEMNodeListView<Dimension>& x) const {
  return std::distance(mDEMNodeListViews.begin(), std::ranges::find(mDEMNodeListViews, x));
}

//------------------------------------------------------------------------------
// calculated appropriates H from a given radius for each DEM nodelist
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
void
DataBaseView<Dimension>::setDEMHfieldFromParticleRadius(const size_t startUniqueIndex,
                                                        const Scalar kernelExtent) {
  for (auto& x: mDEMNodeListViews) x.setHfieldFromParticleRadius(startUniqueIndex,
                                                                 kernelExtent);
}

//------------------------------------------------------------------------------
// The maximum kernel extent being used for DEM
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
typename Dimension::Scalar
DataBaseView<Dimension>::maxDEMNeighborSearchBuffer() const {
  Scalar result = 0.0;
  for (const auto& x: mDEMNodeListViews) result = std::max(result, x.neighborSearchBuffer());
  return result;
  // return std::ranges::max_element(mDEMNodeListViews, {},
  //                                 [](const auto& x) { return x.neighborSearchBuffer(); });
}

}

