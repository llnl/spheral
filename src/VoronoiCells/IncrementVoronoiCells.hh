//---------------------------------Spheral++----------------------------------//
// IncrementVoronoiCells
//
// Specialization of UpdatePolicyBase to advance the Vornoi cell geometry
// during a step without actually recomputing the geometry.  Instead we distort
// the cells by the local velocity gradient.
//
// Created by JMO, Mon May 20 16:04:51 PDT 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_IncrementVoronoiCells_hh__
#define __Spheral_IncrementVoronoiCells_hh__

#include "DataBase/UpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension>
class IncrementVoronoiCells: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using FacetedVolume = typename Dimension::FacetedVolume;

  // Constructors, destructor.
  IncrementVoronoiCells();
  virtual ~IncrementVoronoiCells() {}
  
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
  IncrementVoronoiCells(const IncrementVoronoiCells& rhs) = delete;
  IncrementVoronoiCells& operator=(const IncrementVoronoiCells& rhs) = delete;
};

}

#endif
