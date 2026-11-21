//---------------------------------Spheral++----------------------------------//
// CRKSPHVoidBoundary -- Special boundary class to create void points off of
// the surface.
//
// Created by JMO, Fri Aug 18 10:31:50 PDT 2017
//
// Modified by:
//----------------------------------------------------------------------------//
#ifndef __Spheral_CRKSPHVoidBoundary__
#define __Spheral_CRKSPHVoidBoundary__

#include "Boundary.hh"

namespace Spheral {

template<typename Dimension> class NodeList;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class DataBase;
template<typename Dimension> class Mesh;

template<typename Dimension>
class CRKSPHVoidBoundary : public Boundary<Dimension> {

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
  CRKSPHVoidBoundary(const FieldList<Dimension, int>& surfacePoint,
                     const FieldList<Dimension, std::vector<Vector>>& etaVoidPoints);
  virtual ~CRKSPHVoidBoundary();

  //**********************************************************************
  // All Boundary conditions must provide the following methods:
  virtual void setGhostNodes(NodeList<Dimension>& nodeList) override;

  // For the computed set of ghost nodes, set the positions and H's.
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

  // Find any internal nodes that are in violation of this Boundary.
  virtual void setViolationNodes(NodeList<Dimension>& nodeList) override;

  // For the computed set of nodes in violation of the boundary, bring them
  // back into compliance (for the positions and H's.)
  virtual void updateViolationNodes(NodeList<Dimension>& nodeList) override;
  //**********************************************************************

  // Prevent the Boundary virtual methods from being hidden
  using Boundary<Dimension>::applyGhostBoundary;
  using Boundary<Dimension>::enforceBoundary;

private:
  //--------------------------- Private Interface ---------------------------//
  const FieldList<Dimension, int>& mSurfacePoint;
  const FieldList<Dimension, std::vector<Vector>>& mEtaVoidPoints;
};

}

#endif
