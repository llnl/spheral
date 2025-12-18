//---------------------------------Spheral++----------------------------------//
// PeriodicBoundary -- Apply a Periodic boundary condition to Spheral++
// Fields.
//
// Created by JMO, Wed Apr 19 14:00:50 PDT 2000
//----------------------------------------------------------------------------//

#ifndef PeriodicBoundary_HH
#define PeriodicBoundary_HH

#include "PlanarBoundary.hh"

namespace Spheral {

template<typename Dimension>
class PeriodicBoundary: public PlanarBoundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  // Constructors and destructors.
  PeriodicBoundary();
  PeriodicBoundary(const GeomPlane<Dimension>& plane1,
                   const GeomPlane<Dimension>& plane2);
  virtual ~PeriodicBoundary();

  // Override the default methods for setting ghost nodes.
  virtual void setGhostNodes(NodeList<Dimension>& nodeList) override;
  virtual void updateGhostNodes(NodeList<Dimension>& nodeList) override;

  // Override the default methods for setting violation nodes.
  virtual void setViolationNodes(NodeList<Dimension>& nodeList) override;
  virtual void updateViolationNodes(NodeList<Dimension>& nodeList) override;

  // Override the methods for setting the enter and exit planes.
  virtual const GeomPlane<Dimension>& enterPlane() const override;
  virtual void setEnterPlane(const GeomPlane<Dimension>& enterPlane) override;

  virtual const GeomPlane<Dimension>& exitPlane() const override;
  virtual void setExitPlane(const GeomPlane<Dimension>& exitPlane) override;

  // Override the culling of ghost nodes.
  virtual void cullGhostNodes(const FieldList<Dimension, size_t>& flagSet,
                              FieldList<Dimension, size_t>& old2newIndexMap,
                              std::vector<size_t>& numNodesRemoved) override;

  // Apply the boundary condition to the ghost nodes in the given Field.
  virtual void applyGhostBoundary(FieldBase<Dimension>& fieldBase) const override;
  virtual void applyGhostBoundary(Field<Dimension, FacetedVolume>& field) const override;

  // Enforce the boundary condition on the violation node values in the given Field.
  virtual void enforceBoundary(Field<Dimension, FacetedVolume>& field) const override;

  // Override the base reset method.
  virtual void reset(const DataBase<Dimension>& dataBase) override;

  // Report the number of ghost nodes in this boundary.
  virtual int numGhostNodes() const override;

  virtual std::string label() const override { return "PeriodicBoundary"; }

  // Prevent the Boundary virtual methods from being hidden
  using Boundary<Dimension>::applyGhostBoundary;
  using Boundary<Dimension>::enforceBoundary;

private:
  //--------------------------- Private Interface ---------------------------//
  // The PeriodicBounary actually is made up of two nested class boundary 
  // conditions, which we hide from the user for simplicity of interface.
//   class PeriodicPlanarBoundary;
#include "PeriodicPlanarBoundary.hh"

  PeriodicPlanarBoundary mPlane1Boundary;
  PeriodicPlanarBoundary mPlane2Boundary;
};

}

#endif
