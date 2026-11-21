//---------------------------------Spheral++----------------------------------//
// ThreeOneDimensionBoundary
// A cheezy way to run 1-D simulations using 3-D objects.  This simply ensures
// that Vectors have 0 x & y components, while tensor are diagonal with the
// yy and zz components set equal to the xx value.
//
// Created by JMO, Wed Jul 28 13:38:28 2004
//----------------------------------------------------------------------------//

#ifndef ThreeOneDimensionBoundary_HH
#define ThreeOneDimensionBoundary_HH

#include "Boundary.hh"

namespace Spheral {

template<typename Dimension>
class ThreeOneDimensionBoundary: public Boundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Boundary<Dimension>::BoundaryNodes BoundaryNodes;

  // Constructors and destructors.
  ThreeOneDimensionBoundary();
  virtual ~ThreeOneDimensionBoundary();

  // Select the ghost nodes.
  virtual void setGhostNodes(NodeList<Dimension>& nodeList);
  virtual void updateGhostNodes(NodeList<Dimension>& nodeList);

  // Find the set of nodes in violation of this boundary in the given NodeList.
  // For planar boundaries this is any node that is "behind" the enter plane.
  virtual void setViolationNodes(NodeList<Dimension>& nodeList);
  virtual void updateViolationNodes(NodeList<Dimension>& nodeList);

  // Apply the boundary condition to ghost nodes.
  virtual void applyGhostBoundary(Field<Dimension, Scalar>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, Vector>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, Tensor>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, SymTensor>& field) const;

  // Apply the boundary condition to the violation node values in the given Field.
  virtual void enforceBoundary(Field<Dimension, Scalar>& field) const;
  virtual void enforceBoundary(Field<Dimension, Vector>& field) const;
  virtual void enforceBoundary(Field<Dimension, Tensor>& field) const;
  virtual void enforceBoundary(Field<Dimension, SymTensor>& field) const;

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#endif
