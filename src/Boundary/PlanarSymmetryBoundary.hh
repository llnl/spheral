//---------------------------------Spheral++----------------------------------//
// PlanarSymmetryBoundary
// This boundary condition enforces (2D) planar symmetry on all fields/nodes by 
// creating (slave) ghost nodes parallel to the computational X-Y plane of 
// internal nodes.  This boundary condition may only be used with 3-D 
// simulations.
//
// Created by JNJ.
//----------------------------------------------------------------------------//

#ifndef PlanarSymmetryBoundary_HH
#define PlanarSymmetryBoundary_HH

#include "Boundary.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

class PlanarSymmetryBoundary: public Boundary<Dim<3> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::Tensor Tensor;
  typedef Dim<3>::SymTensor SymTensor;

  typedef Boundary<Dim<3> >::BoundaryNodes BoundaryNodes;

  // Constructors and destructors.
  PlanarSymmetryBoundary(TableKernel<Dim<3> >* kernel);
  virtual ~PlanarSymmetryBoundary();

  // Select the ghost nodes.
  virtual void setGhostNodes(NodeList<Dim<3> >& nodeList);
  virtual void updateGhostNodes(NodeList<Dim<3> >& nodeList);

  // Find the set of nodes in violation of this boundary in the given NodeList.
  // For planar boundaries this is any node that is "behind" the enter plane.
  virtual void setViolationNodes(NodeList<Dim<3> >& nodeList);
  virtual void updateViolationNodes(NodeList<Dim<3> >& nodeList);

  // Apply the boundary condition to ghost nodes.
  virtual void applyGhostBoundary(Field<Dim<3> , Scalar>& field) const;
  virtual void applyGhostBoundary(Field<Dim<3> , Vector>& field) const;
  virtual void applyGhostBoundary(Field<Dim<3> , Tensor>& field) const;
  virtual void applyGhostBoundary(Field<Dim<3> , SymTensor>& field) const;

  // Apply the boundary condition to the violation node values in the given Field.
  virtual void enforceBoundary(Field<Dim<3> , Scalar>& field) const;
  virtual void enforceBoundary(Field<Dim<3> , Vector>& field) const;
  virtual void enforceBoundary(Field<Dim<3> , Tensor>& field) const;
  virtual void enforceBoundary(Field<Dim<3> , SymTensor>& field) const;

private:

  TableKernel<Dim<3> >* mKernel;
  size_t mFirstGhostIndex;
  double mL, mdz;
};

}

#endif
