//---------------------------------Spheral++----------------------------------//
// PythonNeighbor
//
// PythonNeighbor just provides the virtual function hooks to support the
// NeighborView interface.  The virutal Neighbor interface is passed along
// for descendants to provide.
//
// NOTE: PythonNeighbor and its descendants are CPU-only (no GPU capability).
//----------------------------------------------------------------------------//

#include "Neighbor/PythonNeighbor.hh"
#include "DataBase/DataBase.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Construct a non-owning callback adapter.
//------------------------------------------------------------------------------
template<typename Dimension>
PythonNeighborCallView<Dimension>::
PythonNeighborCallView(PythonNeighbor<Dimension>* parent):
  mParent(parent) {
  REQUIRE(parent != nullptr);
}

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
PythonNeighbor<Dimension>::
PythonNeighbor(DataBase<Dimension>& dataBase,
               const NeighborSearchType searchType,
               const Scalar kernelExtent):
  Neighbor<Dimension>(dataBase),
  NeighborView<Dimension, PythonNeighbor<Dimension>>(searchType, kernelExtent) {
}

//------------------------------------------------------------------------------
// Force the nodeExtent FieldList to be sized correctly
//------------------------------------------------------------------------------
template<typename Dimension>
void
PythonNeighbor<Dimension>::
resizeNodeExtents() {
  mDataBase.resizeGlobalFieldList(mNodeExtent,
                                  Vector::zero(),
                                  "node extents",
                                  false);
  mNodeExtentView = mNodeExtent.view();
}

//------------------------------------------------------------------------------
// Clear our state
//------------------------------------------------------------------------------
template<typename Dimension>
void
PythonNeighbor<Dimension>::
clear() {
  resizeNodeExtents();
}

//------------------------------------------------------------------------------
// PythonNeighborCallView forwarding methods.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PythonNeighborCallView<Dimension>::
resizeNodeExtents() {
  mParent->resizeNodeExtents();
}

template<typename Dimension>
void
PythonNeighborCallView<Dimension>::
clear() {
  mParent->clear();
}

template<typename Dimension>
void
PythonNeighborCallView<Dimension>::
updateNodes() {
  mParent->updateNodes();
}

template<typename Dimension>
void
PythonNeighborCallView<Dimension>::
updateNodes(const std::span<NodeID>& nodeIDs) {
  mParent->updateNodes(nodeIDs);
}

template<typename Dimension>
void
PythonNeighborCallView<Dimension>::
reinitialize() {
  mParent->reinitialize();
}

template<typename Dimension>
void
PythonNeighborCallView<Dimension>::
reinitialize(const Vector& xmin,
             const Vector& xmax,
             const Scalar htarget) {
  mParent->reinitialize(xmin, xmax, htarget);
}

template<typename Dimension>
const DataBase<Dimension>&
PythonNeighborCallView<Dimension>::
dataBase() const {
  return mParent->dataBase();
}

template<typename Dimension>
const FieldList<Dimension, typename Dimension::Vector>&
PythonNeighborCallView<Dimension>::
nodeExtent() const {
  return mParent->nodeExtent();
}

template<typename Dimension>
const FieldListView<Dimension, typename Dimension::Vector>&
PythonNeighborCallView<Dimension>::
nodeExtentView() const {
  return mParent->nodeExtentView();
}

template<typename Dimension>
NeighborSearchType
PythonNeighborCallView<Dimension>::
neighborSearchType() const {
  return mParent->neighborSearchType();
}

template<typename Dimension>
typename Dimension::Scalar
PythonNeighborCallView<Dimension>::
kernelExtent() const {
  return mParent->kernelExtent();
}

template<typename Dimension>
void
PythonNeighborCallView<Dimension>::
setNodeExtents(const DataBaseView<Dimension>& dataBase) {
  mParent->setNodeExtents(dataBase);
}

template<typename Dimension>
void
PythonNeighborCallView<Dimension>::
setInternalNodeExtents(const DataBaseView<Dimension>& dataBase) {
  mParent->setInternalNodeExtents(dataBase);
}

template<typename Dimension>
void
PythonNeighborCallView<Dimension>::
setGhostNodeExtents(const DataBaseView<Dimension>& dataBase) {
  mParent->setGhostNodeExtents(dataBase);
}

template<typename Dimension>
void
PythonNeighborCallView<Dimension>::
setNodeExtents(const DataBaseView<Dimension>& dataBase,
               const Container<NodeID>& nodeIDs) {
  mParent->setNodeExtents(dataBase, nodeIDs);
}

template<typename Dimension>
size_t
PythonNeighborCallView<Dimension>::
precullList(const FieldListView<Dimension, Vector>& positions,
            const Vector& minMasterPosition,
            const Vector& maxMasterPosition,
            const Vector& minMasterExtent,
            const Vector& maxMasterExtent,
            const size_t numCoarse,
            Container<NodeID>& coarseList) const {
  return mParent->precullList(positions,
                              minMasterPosition, maxMasterPosition,
                              minMasterExtent, maxMasterExtent,
                              numCoarse, coarseList);
}

template<typename Dimension>
typename Dimension::Vector
PythonNeighborCallView<Dimension>::
Hextent(const Scalar& H, const Scalar& kernelExtent) {
  return PythonNeighbor<Dimension>::Hextent(H, kernelExtent);
}

template<typename Dimension>
typename Dimension::Vector
PythonNeighborCallView<Dimension>::
Hextent(const SymTensor& H, const Scalar& kernelExtent) {
  return PythonNeighbor<Dimension>::Hextent(H, kernelExtent);
}

template<typename Dimension>
std::tuple<bool, size_t, size_t>
PythonNeighborCallView<Dimension>::
setMasterList(const Vector& position,
              const SymTensor& H,
              Container<NodeID>& masterList,
              Container<NodeID>& coarseNeighbors,
              const bool ghostConnectivity) const {
  return mParent->setMasterList(position, H, masterList, coarseNeighbors,
                                ghostConnectivity);
}

template<typename Dimension>
size_t
PythonNeighborCallView<Dimension>::
setRefineNeighborList(const Vector& position,
                      const SymTensor& H,
                      const size_t numCoarse,
                      const Container<NodeID>& coarseNeighbors,
                      Container<NodeID>& refineNeighbors) const {
  return mParent->setRefineNeighborList(position, H, numCoarse,
                                        coarseNeighbors, refineNeighbors);
}

template<typename Dimension>
std::tuple<bool, size_t, size_t>
PythonNeighborCallView<Dimension>::
setMasterList(const Vector& position,
              Container<NodeID>& masterList,
              Container<NodeID>& coarseNeighbors,
              const bool ghostConnectivity) const {
  return mParent->setMasterList(position, masterList, coarseNeighbors,
                                ghostConnectivity);
}

template<typename Dimension>
size_t
PythonNeighborCallView<Dimension>::
setRefineNeighborList(const Vector& position,
                      const size_t numCoarse,
                      const Container<NodeID>& coarseNeighbors,
                      Container<NodeID>& refineNeighbors) const {
  return mParent->setRefineNeighborList(position, numCoarse,
                                        coarseNeighbors, refineNeighbors);
}

template<typename Dimension>
std::tuple<bool, size_t, size_t>
PythonNeighborCallView<Dimension>::
setMasterList(const GeomPlane<Dimension>& enterPlane,
              const GeomPlane<Dimension>& exitPlane,
              Container<NodeID>& masterList,
              Container<NodeID>& coarseNeighbors) const {
  return mParent->setMasterList(enterPlane, exitPlane,
                                masterList, coarseNeighbors);
}

}
