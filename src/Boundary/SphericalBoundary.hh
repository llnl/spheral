//---------------------------------Spheral++----------------------------------//
// SphericalBoundary -- Create a 3-D boundary around a line of nodes down
// the x-axis appropriate for enforcing a spherically symmetric system.
//
// Created by JMO, Tue Mar 15 21:39:43 PST 2005
//----------------------------------------------------------------------------//

#ifndef SphericalBoundary_HH
#define SphericalBoundary_HH

#include "Boundary.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "Geometry/Dimension.hh"
#include "DataOutput/registerWithRestart.hh"

namespace Spheral {

class FileIO;

class SphericalBoundary: public Boundary<Dim<3> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<3> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;
  typedef Dimension::ThirdRankTensor ThirdRankTensor;
  typedef Dimension::FourthRankTensor FourthRankTensor;
  typedef Dimension::FifthRankTensor FifthRankTensor;
  typedef Dimension::FacetedVolume FacetedVolume;
  typedef Boundary<Dimension>::BoundaryNodes BoundaryNodes;

  // Constructors and destructors.
  SphericalBoundary(const DataBase<Dim<3> >& dataBase);
  virtual ~SphericalBoundary();

  // Use the given NodeList's neighbor object to select the ghost nodes.
  virtual void setGhostNodes(NodeList<Dimension>& nodeList) override;
  virtual void updateGhostNodes(NodeList<Dimension>& nodeList) override;

  // Apply the boundary condition to the ghost node values in the given Field.
  virtual void applyGhostBoundary(Field<Dimension, int>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, Scalar>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, Vector>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, Tensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, SymTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, ThirdRankTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, FourthRankTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, FifthRankTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, std::vector<Scalar> >& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, FacetedVolume>& field) const override;

  // Find the set of nodes in violation of this boundary in the given NodeList.
  virtual void setViolationNodes(NodeList<Dimension>& nodeList) override;
  virtual void updateViolationNodes(NodeList<Dimension>& nodeList) override;

  // Apply the boundary condition to the violation node values in the given Field.
  virtual void enforceBoundary(Field<Dimension, int>& field) const override;
  virtual void enforceBoundary(Field<Dimension, Scalar>& field) const override;
  virtual void enforceBoundary(Field<Dimension, Vector>& field) const override;
  virtual void enforceBoundary(Field<Dimension, Tensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, SymTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, ThirdRankTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, FourthRankTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, FifthRankTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, FacetedVolume>& field) const override;

  // Find the effective reflection operator between the given x-axis
  // position and slaved ghost position.
  Tensor reflectOperator(const Vector& r0, const Vector& r1) const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "SphericalBoundary"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

  // Prevent the Boundary virtual methods from being hidden
  using Boundary<Dimension>::applyGhostBoundary;
  using Boundary<Dimension>::enforceBoundary;

private:
  //--------------------------- Private Interface ---------------------------//
  FieldList<Dim<3>, std::vector<Vector> > mGhostPositions;

  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#include "SphericalBoundaryInline.hh"

#endif
