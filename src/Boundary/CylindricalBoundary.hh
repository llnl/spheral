//---------------------------------Spheral++----------------------------------//
// CylindricalBoundary -- Create a 3-D boundary around a set of nodes in the
// xy plane to emulate cylindrical (RZ) coordinates.
//
// Created by JMO, Tue Mar 29 13:48:03 PST 2005
//----------------------------------------------------------------------------//

#ifndef CylindricalBoundary_HH
#define CylindricalBoundary_HH

#include "Boundary.hh"
#include "Geometry/Dimension.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "DataOutput/registerWithRestart.hh"

namespace Spheral {

class FileIO;

class CylindricalBoundary: public Boundary<Dim<3> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Dimension = Dim<3>;
  using Scalar = Dimension::Scalar;
  using Vector = Dimension::Vector;
  using Tensor = Dimension::Tensor;
  using SymTensor = Dimension::SymTensor;
  using ThirdRankTensor = Dimension::ThirdRankTensor;
  using FourthRankTensor = Dimension::FourthRankTensor;
  using FifthRankTensor = Dimension::FifthRankTensor;
  using FacetedVolume = Dimension::FacetedVolume;

  using BoundaryNodes = Boundary<Dimension>::BoundaryNodes;

  // Constructors and destructors.
  CylindricalBoundary(const DataBase<Dim<3> >& dataBase);
  virtual ~CylindricalBoundary() = default;

  // Select the ghost nodes.
  virtual void setGhostNodes(NodeList<Dimension>& nodeList) override;
  virtual void updateGhostNodes(NodeList<Dimension>& nodeList) override;

  // Apply the boundary condition to the ghost node values in the given Field.
  virtual void applyGhostBoundary(Field<Dimension, Vector>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, Tensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, SymTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, ThirdRankTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, FourthRankTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, FifthRankTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, FacetedVolume>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, RKCoefficients<Dimension>>& field) const override;

  // Find the set of nodes in violation of this boundary in the given NodeList.
  virtual void setViolationNodes(NodeList<Dimension>& nodeList) override;
  virtual void updateViolationNodes(NodeList<Dimension>& nodeList) override;

  // Apply the boundary condition to the violation node values in the given Field.
  virtual void enforceBoundary(Field<Dimension, Vector>& field) const override;
  virtual void enforceBoundary(Field<Dimension, Tensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, SymTensor>& field) const override;

  // Find the effective reflection operator between the given xy
  // position and slaved ghost position.
  static Tensor reflectOperator(const Vector& r0, const Vector& r1);

  // Compute the target angular spacing for the given position and H.
  static double angularSpacing(const double ri,
                               const double hzi,
                               const double nodePerh,
                               const double kernelExtent);

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "CylindricalBoundary"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

  // Prevent the Boundary virtual methods from being hidden
  using Boundary<Dimension>::applyGhostBoundary;
  using Boundary<Dimension>::enforceBoundary;

private:
  //--------------------------- Private Interface ---------------------------//
  const DataBase<Dim<3>>& mDataBase;
  FieldList<Dim<3>, Scalar> mDeltaPhi;
  FieldList<Dim<3>, Vector> mGhostPositions;

  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#include "CylindricalBoundaryInline.hh"

#endif
