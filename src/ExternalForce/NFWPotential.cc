//---------------------------------Spheral++----------------------------------//
// NFWPotential -- Impose a potential due to a Navarro Frenk White mass 
// density profile, which is supposed to be a "universal" dark matter halo
// profile seen in cosmological N body simulations.
// Navarro, Frenk, & White 1997, ApJ, 490, 493-508.
//
// Note that this object is templated on your unit choice, since we depend on
// G and time units (via how the Hubble constant is defined.)
//
// NFW density profile:
// rho/rho_0(r) = delta_c/{ (r/r_s)*(1 + r/r_s)^2 }
//
// Corresponding mass profile:
// M(r) = 4/3 pi delta_c rho_0 r_s^4 { 1 + r/r_s - 2 ln(1 + r/r_s) - 1/(1 + r/r_s) }
//
// Balancing velocity profile:
// v(r) = 4/3 pi delta_c rho_0 r_s^4 G {( 1/(1 + r/r_s) + 2 ln(1 + r/r_s) - 1 )/r +
//                                      1/(1 + r/r_s)^2 - 2/(1 + r/r_s) }
//
// Created by JMO, Thu May  8 18:02:30 PDT 2003
//----------------------------------------------------------------------------//
#include "NFWPotential.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

#include <limits>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
NFWPotential<Dimension>::
NFWPotential(double deltac, double rs, double h0, 
             const typename Dimension::Vector& origin,
             const PhysicalConstants& constants):
  GenericBodyForce<Dimension>(),
  mDeltac(deltac),
  mRs(rs),
  mh0(h0),
  mOrigin(origin),
  mConstants(constants),
  mDeltaPhiFraction(0.01),
  mCriticalDensity(0.0),
  mPotentialEnergy(0.0) {
  ENSURE(deltac > 0.0);
  ENSURE(rs > 0.0);
  ENSURE(h0 > 0.0);

  // Set the hubble constant, which also sets the critical density.
  seth0(h0);
  CHECK(mCriticalDensity > 0.0);
}

//------------------------------------------------------------------------------
// Calculate the acceleration due to the potential of the profile.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NFWPotential<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Get the node positions from the state.
  const auto mnode = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  const auto numNodeLists = mnode.numFields();
  CHECK(position.numFields() == numNodeLists);
  CHECK(velocity.numFields() == numNodeLists);

  // Get the acceleration and position change vectors we'll be modifying.
  auto DxDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero());
  auto DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero());

  // Loop over the internal nodes.
  mPotentialEnergy = 0.0;
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto n = mnode[k]->numInternalElements();
    for (auto i = 0u; i < n; ++i) {
      const auto& r = position(k, i) - mOrigin;
      const auto  runit = r.unitVector();
      const auto  rsoft2 = r.magnitude2() + 1.0e-10;
      const auto  rsoft = sqrt(rsoft2);
      const auto  Mr = enclosedMass(rsoft);
      const auto  thpt = mConstants.G()*Mr;
      DxDt(k, i) += velocity(k, i);
      DvDt(k, i) -= thpt/rsoft2*runit;
      mPotentialEnergy -= thpt*mnode(k, i)/rsoft;
    }
  }
}

//------------------------------------------------------------------------------
// Calculate the timestep constraint.
//------------------------------------------------------------------------------
template<typename Dimension>
typename NFWPotential<Dimension>::TimeStepType
NFWPotential<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& /*derivs*/,
   const typename Dimension::Scalar /*currentTime*/) const {

  // This is a hard one to pick for this package.  We'll choose a timestep
  // such that no nodes potential should change more than a set fraction.
  // (Delta phi)/phi = -(Delta r)/(r^2 + rc^2)
  // dt = (Delta phi)/phi (r^2 + rc^2) / v
  Scalar mindt = std::numeric_limits<Scalar>::max();
  Scalar minr, minv;
  const auto position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  const auto numNodeLists = position.numFields();
  CHECK(velocity.numFields() == numNodeLists);
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto n = position[k]->numInternalElements();
    for (auto i = 0u; i < n; ++i) {
      const auto r = position(k, i) - mOrigin;
      const auto v = velocity(k, i).magnitude() + 1.0e-10;
      const auto rsoft2 = r.magnitude2() + 1.0e-10;
      const auto rsoft = sqrt(rsoft2);
      const auto rsoft3 = rsoft*rsoft2;
      const auto r0 = rsoft/mRs;
      const auto r1 = 1.0 + r0;
      const auto phi = (1.0/r1 + 2.0*log(r1) - 1.0)/rsoft2 + 1.0/(r1*r1*rsoft) - 2.0/(r1*rsoft);
      const auto dPhidr = -r1/(rsoft2*r1*r1*mRs) - 2.0/(rsoft3*r1) +
        2.0/(rsoft2*r1*mRs) - 4.0*log(r1)/rsoft3 +
        2.0/rsoft3 - 2.0/(r1*r1*r1*mRs*rsoft) -
        1.0/(rsoft2*r1*r1) + 2.0/(rsoft*r1*r1*mRs) +
        2.0/(rsoft2*r1);
      const auto dt = mDeltaPhiFraction*phi/(dPhidr*v);
      if (dt < mindt) {
        mindt = dt;
        minr = rsoft;
        minv = v;
      }
    }
  }

  std::stringstream reasonStream;
  reasonStream << "mindt = " << mindt << " | "
               << "rsoft = " << minr << " " 
               << "minv = " << minv << std::endl;

  return TimeStepType(mindt, reasonStream.str());
}

}
