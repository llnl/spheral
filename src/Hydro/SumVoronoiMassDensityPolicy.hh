//---------------------------------Spheral++----------------------------------//
// SumVoronoiMassDensityPolicy -- An implementation of UpdatePolicyBase 
// specialized for the updating the mass density according to the specific 
// volume from the Voronoi tesselation.  This version uses a weighted sum of the
// local mass and volume to do the deed.
//
// Created by JMO, Mon Aug  1 10:48:03 PDT 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_SumVoronoiMassDensityPolicy_hh__
#define __Spheral_SumVoronoiMassDensityPolicy_hh__

#include <string>

#include "DataBase/UpdatePolicyBase.hh"
#include "Kernel/TableKernel.hh"
#include "Physics/Physics.hh"

namespace Spheral {

template<typename Dimension>
class SumVoronoiMassDensityPolicy: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  SumVoronoiMassDensityPolicy(const TableKernel<Dimension>& W,
                              const Physics<Dimension>& package,
                              const double rhoMin,
                              const double rhoMax);
  virtual ~SumVoronoiMassDensityPolicy() = default;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Forbidden methods
  SumVoronoiMassDensityPolicy(const SumVoronoiMassDensityPolicy& rhs) = delete;
  SumVoronoiMassDensityPolicy& operator=(const SumVoronoiMassDensityPolicy& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  bool mEnforceBoundaries;
  const TableKernel<Dimension>& mW;
  const Physics<Dimension>& mPackage;
  double mRhoMin, mRhoMax;
};

}

#endif
