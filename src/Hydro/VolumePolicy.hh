//---------------------------------Spheral++----------------------------------//
// VolumePolicy -- An implementation of ReplaceState specialized
// for the updating the volume based on the current Voronoi tesselation.
//
// Created by JMO, Mon Aug  1 15:17:36 PDT 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_VolumePolicy_hh__
#define __Spheral_VolumePolicy_hh__

#include <string>

#include "DataBase/UpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension>
class VolumePolicy: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  VolumePolicy();
  virtual ~VolumePolicy() = default;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // We'll make the updateAsIncrement a no-op.
  virtual void updateAsIncrement(const KeyType& /*key*/,
                                 State<Dimension>& /*state*/,
                                 StateDerivatives<Dimension>& /*derivs*/,
                                 const double /*multiplier*/,
                                 const double /*t*/,
                                 const double /*dt*/) override {}

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Forbidden methods
  VolumePolicy(const VolumePolicy& rhs) = delete;
  VolumePolicy& operator=(const VolumePolicy& rhs) = delete;
};

}

#endif
