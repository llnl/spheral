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
// v^2(r) = 4/3 pi delta_c rho_0 r_s^4 G {( 1/(1 + r/r_s) + 2 ln(1 + r/r_s) - 1 )/r +
//                                        1/(1 + r/r_s)^2 - 2/(1 + r/r_s) }
//
// Created by JMO, Thu May  8 18:02:30 PDT 2003
//----------------------------------------------------------------------------//
#ifndef NFWPotential_HH
#define NFWPotential_HH

#include "Physics/GenericBodyForce.hh"
#include "Material/PhysicalConstants.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class NodeList;
template<typename Dimension> class DataBase;
class PhysicalConstants;

template<typename Dimension>
class NFWPotential: public GenericBodyForce<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using TimeStepType = typename Physics<Dimension>::TimeStepType;

  // Constructors.
  // deltac = characteristic dimensionless density
  // rs = scale radius
  // h0 = Hubble constant (units of 100 km/sec/Mpc)
  NFWPotential(double deltac,
               double rs,
               double h0,
               const Vector& origin,
               const PhysicalConstants& constants);

  // Destructor.
  virtual ~NFWPotential() = default;

  // This is the derivative method that all BodyPotential classes must provide.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const override;

  // Provide the timestep appropriate for this package.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  // Get the cumulative potential energy calculated in the last 
  // evaluateDerivatives.
  virtual Scalar extraEnergy() const override;

  //! Required label for Physics interface.
  virtual std::string label() const override { return "NFWPotential"; }

  // Calculate the mass density as a function of radius.
  Scalar massDensity(Scalar r) const;

  // Calculate the enclosed mass as a function of radius.
  Scalar enclosedMass(Scalar r) const;

  // Calculate the supporting orbital velocity as a function of radius.
  Scalar orbitalVelocity(Scalar r) const;

  // The characteristic density.
  Scalar characteristicDensity() const;
  void setCharacteristicDensity(const Scalar x);

  // The scale Radius.
  Scalar scaleRadius() const;
  void setScaleRadius(const Scalar x);

  // The hubble constant (in units of 100 km/sec/Mpc).
  Scalar h0() const;
  void seth0(const Scalar rc);

  // The critical (closure) density for the universe.
  Scalar criticalDensity() const;

  // The origin (center of the profile).
  const Vector& origin() const;
  void setOrigin(const Vector& origin);

  // The maximum allowed fractional change in a particles potential (for 
  // setting the timestep).
  Scalar deltaPotentialFraction() const;
  void setDeltaPotentialFraction(const Scalar deltaPhi);

  // No default constructor, copying, or assignment.
  NFWPotential() = delete;
  NFWPotential(const NFWPotential& rhs) = delete;
  NFWPotential& operator=(const NFWPotential& rhs) = delete;

private:
  //--------------------------- Public Interface ---------------------------//
  Scalar mDeltac;
  Scalar mRs;
  Scalar mh0;
  Vector mOrigin;
  PhysicalConstants mConstants;
  Scalar mDeltaPhiFraction;
  Scalar mCriticalDensity;
  mutable Scalar mPotentialEnergy;
};

}

#include "NFWPotentialInline.hh"

#endif
