//---------------------------------Spheral++----------------------------------//
// UpdateVoronoiCells
//
// Specialization of UpdatePolicyBase to advance the Vornoi cell geometry
// during a step without actually recomputing the geometry.  Instead we distort
// the cells by the local velocity gradient.
//
// Created by JMO, Mon May 20 16:04:51 PDT 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_UpdateVoronoiCells_hh__
#define __Spheral_UpdateVoronoiCells_hh__

#include "DataBase/UpdatePolicyBase.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"

namespace Spheral {

template<typename Dimension>
class UpdateVoronoiCells: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using FacetedVolume = typename Dimension::FacetedVolume;

  // Constructors, destructor.
  UpdateVoronoiCells(FieldList<Dimension, Scalar>& volume,
                     FieldList<Dimension, Scalar>& weight,
                     FieldList<Dimension, Vector>& deltaCentroid,
                     FieldList<Dimension, vector<Vector>>& etaVoidPoints,
                     const std::vector<Boundary<Dimension>*>& boundaries,
                     const std::vector<FacetedVolume>& facetedBoundaries,
                     const std::vector<std::vector<FacetedVolume>>& facetedHoles);
  virtual ~UpdateVoronoiCells() {}
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // No default constructor or copying
  UpdateVoronoiCells(const UpdateVoronoiCells& rhs) = delete;
  UpdateVoronoiCells& operator=(const UpdateVoronoiCells& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  FieldList<Dimension, Scalar>& mVolume;
  FieldList<Dimension, Scalar>& mWeight;
  FieldList<Dimension, Vector>& mDeltaCentroid;
  FieldList<Dimension, std::vector<Vector>>& mEtaVoidPoints;
  const std::vector<Boundary<Dimension>*>& mBoundaries;
  const std::vector<FacetedVolume>& mFacetedBoundaries;
  const std::vector<std::vector<FacetedVolume>>& mFacetedHoles;
};

}

#endif
