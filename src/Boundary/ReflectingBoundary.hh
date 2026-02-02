//---------------------------------Spheral++----------------------------------//
// ReflectingBoundary -- Apply a Reflecting boundary condition to Spheral++
// Fields.
//
// Created by JMO, Wed Feb 16 21:01:02 PST 2000
//----------------------------------------------------------------------------//

#ifndef ReflectingBoundary_HH
#define ReflectingBoundary_HH

#include "Boundary/Boundary.hh"
#include "Boundary/PlanarBoundary.hh"
#include "RK/RKCorrectionParams.hh"

#include "Eigen/Sparse"

#include <unordered_map>

namespace Spheral {

template<typename Dimension>
class ReflectingBoundary: public PlanarBoundary<Dimension> {

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
  typedef typename Eigen::SparseMatrix<double> TransformationMatrix;

  // Constructors and destructors.
  ReflectingBoundary();
  ReflectingBoundary(const GeomPlane<Dimension>& plane);
  virtual ~ReflectingBoundary();

  // Apply the boundary condition to the ghost values of given Field.
  virtual void applyGhostBoundary(Field<Dimension, Vector>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, Tensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, SymTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, ThirdRankTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, FourthRankTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, FifthRankTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, FacetedVolume>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, RKCoefficients<Dimension>>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, std::vector<Vector>>& field) const override;

  // Apply the boundary condition to the violation node values in the given Field.
  virtual void enforceBoundary(Field<Dimension, Vector>& field) const override;
  virtual void enforceBoundary(Field<Dimension, Tensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, SymTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, ThirdRankTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, FourthRankTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, FifthRankTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, FacetedVolume>& field) const override;
  virtual void enforceBoundary(Field<Dimension, RKCoefficients<Dimension>>& field) const override;
  virtual void enforceBoundary(Field<Dimension, std::vector<Vector>>& field) const override;

  // Apply the boundary condition to face centered fields on a tessellation.
  virtual void enforceBoundary(std::vector<int>& faceField, const Mesh<Dimension>& mesh) const override;
  virtual void enforceBoundary(std::vector<Scalar>& faceField, const Mesh<Dimension>& mesh) const override;
  virtual void enforceBoundary(std::vector<Vector>& faceField, const Mesh<Dimension>& mesh) const override;
  virtual void enforceBoundary(std::vector<Tensor>& faceField, const Mesh<Dimension>& mesh) const override;
  virtual void enforceBoundary(std::vector<SymTensor>& faceField, const Mesh<Dimension>& mesh) const override;
  virtual void enforceBoundary(std::vector<ThirdRankTensor>& faceField, const Mesh<Dimension>& mesh) const override;
  virtual void enforceBoundary(std::vector<FourthRankTensor>& faceField, const Mesh<Dimension>& mesh) const override;
  virtual void enforceBoundary(std::vector<FifthRankTensor>& faceField, const Mesh<Dimension>& mesh) const override;

  // Fill in faces on this boundary with effective opposite face values.
  virtual void swapFaceValues(Field<Dimension, std::vector<Scalar> >& /*field*/,
                              const Mesh<Dimension>& /*mesh*/) const override;
  virtual void swapFaceValues(Field<Dimension, std::vector<Vector> >& field,
                              const Mesh<Dimension>& mesh) const override;

  // Allow read only access to the reflection operators.
  const Tensor& reflectOperator() const;
  const TransformationMatrix& rkReflectOperator(const RKOrder order,
                                                const bool useHessian) const;

  // Valid test.
  virtual bool valid() const override;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "ReflectingBoundary"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //****************************************************************************

  // Prevent the Boundary virtual methods from being hidden
  using Boundary<Dimension>::applyGhostBoundary;
  using Boundary<Dimension>::enforceBoundary;

private:
  //--------------------------- Private Interface ---------------------------//
  Tensor mReflectOperator;
  std::unordered_map<RKOrder, std::pair<TransformationMatrix, TransformationMatrix>> mrkReflectOperators;
};

}

#include "ReflectingBoundaryInline.hh"

#endif
