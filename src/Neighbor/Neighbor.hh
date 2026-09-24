//---------------------------------Spheral++----------------------------------//
// Neighbor -- An abstract base class for the Neighbor objects.
//
// Created by J. Michael Owen, Sun Nov 12 10:33:55 2000
// Split into Value/View pattern Tue Sep 22 13:55:08 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_Neighbor__
#define __Spheral_Neighbor__

#include "Neighbor/NeighborBase.hh"
#include "Geometry/GeomPlane.hh"
#include "Field/Field.hh"
#include "chai/managed_ptr.hpp"

namespace Spheral {

template<typename Dimension> class NodeList;
template<typename Dimension> class NeighborView;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class Neighbor: public NeighborBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using ViewType = NeighborView<Dimension>;

  // Constructors and destructors
  Neighbor(NodeList<Dimension>& nodeList, 
           const NeighborSearchType searchType,
           const Scalar kernelExtent);
  virtual ~Neighbor();

  // Access to the field of node extents.
  const Field<Dimension, Vector>& nodeExtentField() const      { return mNodeExtent; }

  // Manage the node list.
  const NodeList<Dimension>& nodeList() const                  { CHECK(mNodeListPtr); return *mNodeListPtr; }
  const NodeList<Dimension>* nodeListPtr() const               { return mNodeListPtr; }
  void nodeList(NodeList<Dimension>& nodeList);
  void unregisterNodeList();

  // Force the update of internal data for the NodeList.
  virtual void updateNodes() = 0;
  virtual void updateNodes(const std::vector<int>& nodeIDs) = 0;

  //...........................................................................
  // Optional hooks to reinitialize based on the desired target h and box dimensions.
  virtual void reinitialize() {};
  virtual void reinitialize(const Vector& xmin,
                            const Vector& xmax,
                            const Scalar& htarget) {};

  //...........................................................................
  // Methods for accessing the view (required for descendant)
  virtual chai::managed_ptr<ViewType> getView()  = 0;

protected:
  //-------------------------- Protected Interface --------------------------//
  // Provide read/write access to the node index vectors for descendent classes.
  Field<Dimension, Vector>& accessNodeExtentField()            { return mNodeExtent; }

private:
  //--------------------------- Private Interface ---------------------------//
  using NeighborBase<Dimension>::mSearchType;
  using NeighborBase<Dimension>::mKernelExtent;
  NodeList<Dimension>* mNodeListPtr;
  Field<Dimension, Vector> mNodeExtent;
};

}

#endif
