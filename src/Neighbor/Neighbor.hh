//---------------------------------Spheral++----------------------------------//
// Neighbor -- An abstract base class for the Neighbor objects.
//
// Created by J. Michael Owen, Sun Nov 12 10:33:55 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_Neighbor__
#define __Spheral_Neighbor__

#include "Neighbor/NeighborBase.hh"
#include "Geometry/GeomPlane.hh"
#include "Field/Field.hh"

#include <vector>

namespace Spheral {

template<typename Dimension> class NodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class Neighbor: public NeighborBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using iterator = std::vector<int>::iterator;
  using const_iterator = std::vector<int>::const_iterator;

  // Constructors and destructors
  Neighbor(NodeList<Dimension>& nodeList, 
           const NeighborSearchType searchType,
           const double kernelExtent);
  virtual ~Neighbor();

  // Choose the type of neighbor search we wish to use.
  NeighborSearchType neighborSearchType() const;
  void neighborSearchType(NeighborSearchType searchType);

  // All neighboring classes need to now how far to sample.
  double kernelExtent() const;
  void kernelExtent(double kernelExtent);

  // Allow access to the field of node extents.
  const Field<Dimension, Vector>& nodeExtentField() const;

  // Access the node list.
  const NodeList<Dimension>& nodeList() const;
  void nodeList(NodeList<Dimension>& nodeList);

  const NodeList<Dimension>* nodeListPtr() const;
  void nodeListPtr(NodeList<Dimension>* nodeListPtr);

  void unregisterNodeList();

  // Determine the node extent for an individual node.
  Vector nodeExtent(int nodeID) const;

  // Force the node extent field to be computed.
  void setNodeExtents();
  void setNodeExtents(const std::vector<int>& nodeIDs);
  void setInternalNodeExtents();
  void setGhostNodeExtents();

  // Set or refine the neighbor lists for a given node ID.
  virtual void setMasterList(int nodeID,
                             std::vector<int>& masterList,
                             std::vector<int>& coarseNeighbors,
                             const bool ghostConnectivity = false) const = 0;
  virtual void setRefineNeighborList(int nodeID,
                                     const std::vector<int>& coarseNeighbors,
                                     std::vector<int>& refineNeighbors) const;

  // Helper method to cull lists of neighbors based on min/max positions and 
  // min/max extents.
  std::vector<int> 
  precullList(const Vector& minMasterPosition, const Vector& maxMasterPostion,
              const Vector& minMasterExtent, const Vector& maxMasterExtent,
              const std::vector<int>& coarseNeighbors) const;

  // ********** Descendent Neighbor types must provide these methods. **********
  // Set or refine the neighbor lists for the given position and smoothing 
  // scale.
  virtual void setMasterList(const Vector& position,
                             const Scalar& H,
                             std::vector<int>& masterList,
                             std::vector<int>& coarseNeighbors,
                             const bool ghostConnectivity = false) const = 0;
  virtual void setMasterList(const Vector& position,
                             const SymTensor& H,
                             std::vector<int>& masterList,
                             std::vector<int>& coarseNeighbors,
                             const bool ghostConnectivity = false) const = 0;

  virtual void setRefineNeighborList(const Vector& position,
                                     const Scalar& H,
                                     const std::vector<int>& coarseNeighbors,
                                     std::vector<int>& refineNeighbors) const = 0;
  virtual void setRefineNeighborList(const Vector& position,
                                     const SymTensor& H,
                                     const std::vector<int>& coarseNeighbors,
                                     std::vector<int>& refineNeighbors) const = 0;

  // Set Neighbors for the given position.
  virtual void setMasterList(const Vector& position,
                             std::vector<int>& masterList,
                             std::vector<int>& coarseNeighbors,
                             const bool ghostConnectivity = false) const = 0;
  virtual void setRefineNeighborList(const Vector& position,
                                     const std::vector<int>& coarseNeighbors,
                                     std::vector<int>& refineNeighbors) const = 0;

  // Set the neighbor lists based on proximity to planes.
  virtual void setMasterList(const GeomPlane<Dimension>& enterPlane,
                             const GeomPlane<Dimension>& exitPlane,
                             std::vector<int>& masterList,
                             std::vector<int>& coarseNeighbors) const = 0;

  // Force the update of internal data for the NodeList.
  virtual void updateNodes() = 0;
  virtual void updateNodes(const std::vector<int>& nodeIDs) = 0;
  //****************************************************************************

  // Optional hooks to reinitialize based on the desired target h and box dimensions.
  virtual void reinitialize() {};
  virtual void reinitialize(const Vector& /*xmin*/, const Vector& /*xmax*/, const Scalar /*htarget*/) {};

  // Determine if the Neighbor is in a valid, ready to use state.
  virtual bool valid() const;

  // A static method for set master/coarse neighbor sets on multiple NodeLists
  // with culling of the coarse set.  This method is static because it is 
  // intended to deal with multiple Neighbor objects (not just this particular
  // one), and must access private data of the Neighbors.
  template<typename NodeListIteratorType>
  static void setMasterNeighborGroup(const Vector& position,
                                     const SymTensor& H,
                                     const NodeListIteratorType& nodeListBegin,
                                     const NodeListIteratorType& nodeListEnd,
                                     const double kernelExtent,
                                     std::vector<std::vector<int>>& masterLists,
                                     std::vector<std::vector<int>>& coarseNeighbors,
                                     const bool ghostConnectivity = false);

  // Determine the maximum extent of a given H smoothing scale along the
  // Cartesian axes.
  static Vector HExtent(const Scalar& H, const double kernelExtent);
  static Vector HExtent(const SymTensor& H, const double kernelExtent);

protected:
  //-------------------------- Protected Interface --------------------------//
  // Provide read/write access to the node index vectors for descendent classes.
  Field<Dimension, Vector>& accessNodeExtentField();

private:
  //--------------------------- Private Interface ---------------------------//
  NeighborSearchType mSearchType;
  double mKernelExtent;
  NodeList<Dimension>* mNodeListPtr;
  Field<Dimension, Vector> mNodeExtent;
};

// We explicitly specialize the HExtent method for 1, 2, & 3 dimensions.
template<> Dim<1>::Vector Neighbor< Dim<1> >::HExtent(const Dim<1>::SymTensor&, const double kernelExtent);
template<> Dim<2>::Vector Neighbor< Dim<2> >::HExtent(const Dim<2>::SymTensor&, const double kernelExtent);
template<> Dim<3>::Vector Neighbor< Dim<3> >::HExtent(const Dim<3>::SymTensor&, const double kernelExtent);

}

#include "NeighborInline.hh"

#endif
