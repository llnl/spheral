//---------------------------------Spheral++----------------------------------//
// PythonNeighborCallView -- A small CPU-only adapter for PythonNeighbor.
//
// The PythonNeighbor is owned elsewhere.  This type is copyable, so it can be
// stored directly in NeighborVariant and NeighborViewVariant without putting
// ownership or a polymorphic object in either variant.
//----------------------------------------------------------------------------//
#ifndef __Spheral_PythonNeighborCallView__
#define __Spheral_PythonNeighborCallView__

#include "Neighbor/NeighborSearchType.hh"
#include "Neighbor/NeighborViewContainer.hh"
#include "Field/FieldListView.hh"

#include <cstddef>
#include <span>
#include <tuple>
#include <utility>

namespace Spheral {

template<typename Dimension> class DataBaseView;
template<typename Dimension> class DataBase;
template<typename Dimension> class GeomPlane;
template<typename Dimension> class PythonNeighbor;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
class PythonNeighborCallView {
public:
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using SymTensor = typename Dimension::SymTensor;
  using NodeID = std::pair<size_t, size_t>;
  template<typename T> using Container = NeighborViewContainer<T>;

  explicit PythonNeighborCallView(PythonNeighbor<Dimension>* parent);

  // Host-side Neighbor value operations.
  void resizeNodeExtents();
  void clear();
  void updateNodes();
  void updateNodes(const std::span<NodeID>& nodeIDs);
  void reinitialize();
  void reinitialize(const Vector& xmin,
                    const Vector& xmax,
                    const Scalar htarget);
  const DataBase<Dimension>& dataBase() const;
  const FieldList<Dimension, Vector>& nodeExtent() const;

  // Common Neighbor view state and operations.
  const FieldListView<Dimension, Vector>& nodeExtentView() const;
  NeighborSearchType neighborSearchType() const;
  Scalar kernelExtent() const;
  void kernelExtent(const Scalar value);
  void setNodeExtents(const DataBaseView<Dimension>& dataBase);
  void setInternalNodeExtents(const DataBaseView<Dimension>& dataBase);
  void setGhostNodeExtents(const DataBaseView<Dimension>& dataBase);
  void setNodeExtents(const DataBaseView<Dimension>& dataBase,
                      const Container<NodeID>& nodeIDs);
  size_t precullList(const FieldListView<Dimension, Vector>& positions,
                     const Vector& minMasterPosition,
                     const Vector& maxMasterPosition,
                     const Vector& minMasterExtent,
                     const Vector& maxMasterExtent,
                     const size_t numCoarse,
                     Container<NodeID>& coarseList) const;
  static Vector Hextent(const Scalar& H, const Scalar& kernelExtent);
  static Vector Hextent(const SymTensor& H, const Scalar& kernelExtent);

  std::tuple<bool, size_t, size_t> setMasterList(const Vector& position,
                                                 const SymTensor& H,
                                                 Container<NodeID>& masterList,
                                                 Container<NodeID>& coarseNeighbors,
                                                 const bool ghostConnectivity = false) const;
  size_t setRefineNeighborList(const Vector& position,
                               const SymTensor& H,
                               const size_t numCoarse,
                               const Container<NodeID>& coarseNeighbors,
                               Container<NodeID>& refineNeighbors) const;
  std::tuple<bool, size_t, size_t> setMasterList(const Vector& position,
                                                 Container<NodeID>& masterList,
                                                 Container<NodeID>& coarseNeighbors,
                                                 const bool ghostConnectivity = false) const;
  size_t setRefineNeighborList(const Vector& position,
                               const size_t numCoarse,
                               const Container<NodeID>& coarseNeighbors,
                               Container<NodeID>& refineNeighbors) const;
  std::tuple<bool, size_t, size_t> setMasterList(const GeomPlane<Dimension>& enterPlane,
                                                 const GeomPlane<Dimension>& exitPlane,
                                                 Container<NodeID>& masterList,
                                                 Container<NodeID>& coarseNeighbors) const;

  PythonNeighborCallView view() const { return *this; }

private:
  PythonNeighbor<Dimension>* mParent;
};

}

#endif
