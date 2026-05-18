//---------------------------------Spheral++----------------------------------//
// Neighbor -- An abstract base class for the Neighbor objects.
//
// Created by J. Michael Owen, Sun Nov 12 10:33:55 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_Neighbor_hh__
#define __Spheral_Neighbor_hh__

#include "Field/Field.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"

#include <stdint.h>
#include <vector>

namespace Spheral {

template<typename Dimension> class NodeList;
template<typename Dimension, typename DataType> class Field;

enum class NeighborSearchType {
  //None = 0,
  Gather = 1,
  Scatter = 2,
  GatherScatter = 3
};

enum NeighborGroupKind {
  FallbackNeighborGroupKind = 0,
  TreeNeighborGroupKind = 1,
  NestedNeighborGroupKind = 2
};

struct NeighborGroupDescriptor {
  int kind = 0;
  int level = 0;
  uint64_t key = 0;

  bool operator==(const NeighborGroupDescriptor& rhs) const {
    return kind == rhs.kind and level == rhs.level and key == rhs.key;
  }

  bool operator<(const NeighborGroupDescriptor& rhs) const {
    return (kind < rhs.kind or
            (kind == rhs.kind and
             (level < rhs.level or
              (level == rhs.level and key < rhs.key))));
  }
};

template<typename Dimension>
class Neighbor {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef std::vector<int>::iterator iterator;
  typedef std::vector<int>::const_iterator const_iterator;

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

  // Count/fill interfaces for preprocessing flat connectivity structures.
  virtual size_t countMasterList(const Vector& position,
                                 const SymTensor& H,
                                 const bool ghostConnectivity = false) const = 0;
  virtual void fillMasterList(const Vector& position,
                              const SymTensor& H,
                              int* result,
                              const bool ghostConnectivity = false) const = 0;
  virtual size_t countCoarseNeighbors(const Vector& position,
                                      const SymTensor& H,
                                      const bool ghostConnectivity = false) const = 0;
  virtual void fillCoarseNeighbors(const Vector& position,
                                   const SymTensor& H,
                                   int* result,
                                   const bool ghostConnectivity = false) const = 0;
  virtual int groupKind() const;
  virtual NeighborGroupDescriptor groupDescriptor(const Vector& position,
                                                  const SymTensor& H,
                                                  const uint64_t fallbackToken) const;

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
