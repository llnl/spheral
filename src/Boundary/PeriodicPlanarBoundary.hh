//---------------------------------Spheral++----------------------------------//
// PeriodicPlanarBoundary -- The nested class member of Periodic Boundary that
// does the actual work.
//
// Created by JMO, Wed Apr 19 15:00:50 PDT 2000
//----------------------------------------------------------------------------//

#ifndef __Spheral_PeriodicPlanarBoundary__
#define __Spheral_PeriodicPlanarBoundary__

class PeriodicPlanarBoundary: public PlanarBoundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ThirdRankTensor = typename Dimension::ThirdRankTensor;
  using FourthRankTensor = typename Dimension::FourthRankTensor;
  using FifthRankTensor = typename Dimension::FifthRankTensor;
  using FacetedVolume = typename Dimension::FacetedVolume;

  // Constructors and destructors.
  PeriodicPlanarBoundary(const DataBase<Dimension>& dataBase,
                         const GeomPlane<Dimension>& plane1,
                         const GeomPlane<Dimension>& plane2);
  virtual ~PeriodicPlanarBoundary() = default;

  // Apply the boundary condition to the ghost nodes in the given Field.
  virtual void applyGhostBoundary(FieldBase<Dimension>& fieldBase) const override;
  virtual void applyGhostBoundary(Field<Dimension, FacetedVolume>& field) const override;

  // Enforce the boundary condition on the violation nodes in the given Field.
  virtual void enforceBoundary(Field<Dimension, FacetedVolume>& field) const override;

  // Valid test.
  virtual bool valid() const override;

private:
  //--------------------------- Private Interface ---------------------------//
};

#endif
