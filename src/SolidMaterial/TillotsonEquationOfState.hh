//---------------------------------Spheral++----------------------------------//
// TillotsonEquationOfState -- This equation of state is designed to represent
// metallic materials over a range of pressure and density -- spanning solid, 
// liquid, and vapor states.
//
// Tillotson 1962
//
// Created by JMO, Wed Mar 16 23:31:17 PDT 2011
//----------------------------------------------------------------------------//
#ifndef TillotsonEquationOfState_HH
#define TillotsonEquationOfState_HH

#include "SolidEquationOfState.hh"

#include <tuple>

namespace Spheral {

template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class TillotsonEquationOfState: public SolidEquationOfState<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  TillotsonEquationOfState(const double referenceDensity,
                           const double etamin,
                           const double etamax,
                           const double etamin_solid,
                           const double etamax_solid,
                           const double a,
                           const double b,
                           const double A,
                           const double B,
                           const double alpha,
                           const double beta,
                           const double eps0,
                           const double epsLiquid,
                           const double epsVapor,
                           const double atomicWeight,
                           const PhysicalConstants& constants,
                           const double externalPressure,
                           const double minimumPressure,
                           const double maximumPressure,
                           const double minimumPressureDamage,
                           const MaterialPressureMinType minPressureType);
  ~TillotsonEquationOfState();

  // We require any equation of state to define the following properties.
  virtual void setPressure(Field<Dimension, Scalar>& Pressure,
                           const Field<Dimension, Scalar>& massDensity,
                           const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setPressureAndDerivs(Field<Dimension, Scalar>& Pressure,           // set pressure
                                    Field<Dimension, Scalar>& dPdu,               // set (\partial P)/(\partial u) (specific thermal energy)
                                    Field<Dimension, Scalar>& dPdrho,             // set (\partial P)/(\partial rho) (density)
                                    const Field<Dimension, Scalar>& massDensity,
                                    const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setTemperature(Field<Dimension, Scalar>& temperature,
                              const Field<Dimension, Scalar>& massDensity,
                              const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                                        const Field<Dimension, Scalar>& massDensity,
                                        const Field<Dimension, Scalar>& temperature) const override;

  virtual void setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                               const Field<Dimension, Scalar>& massDensity,
                               const Field<Dimension, Scalar>& temperature) const override;

  virtual void setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setGammaField(Field<Dimension, Scalar>& gamma,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setEntropy(Field<Dimension, Scalar>& entropy,
                          const Field<Dimension, Scalar>& massDensity,
                          const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  // Access the member data.
  double etamin_solid() const;
  double etamax_solid() const;
  double a() const;
  double b() const;
  double A() const;
  double B() const;
  double alpha() const;
  double beta() const;
  double eps0() const;
  double epsLiquid() const;
  double epsVapor() const;
  double atomicWeight() const;
  
  void etamin_solid(double x);
  void etamax_solid(double x);
  void a(double x);
  void b(double x);
  void A(double x);
  void B(double x);
  void alpha(double x);
  void beta(double x);
  void eps0(double x);
  void epsLiquid(double x);
  void epsVapor(double x);
  void atomicWeight(double x);
  
  // Compute the derivative of the pressure with respect to the density.
  double computeDPDrho(const Scalar massDensity,
                       const Scalar specificThermalEnergy) const;

  // Helper methods for computing the components of the pressure and it's derivatives.
  double computePhi(const double& eta, const double& eps) const;
  double computeP1(const double& mu, const double& P2) const;
  double computeP2(const double& phi, const double& mu, const double& rho, const double& eps) const;
  double computeP4(const double& phi, const double& mu, const double& eta, const double& rho, const double& eps) const;
  double compute_dphidrho_eps(const double& rho0, const double& eta, const double& eps) const;
  double compute_dP1drho_eps(const double& rho0, const double& mu, const double& dP2drho_eps) const;
  double compute_dP2drho_eps(const double& phi, const double& dphidrho_eps, const double& rho0, const double& rho, const double& eps) const;
  double compute_dP4drho_eps(const double& phi, const double& dphidrho_eps, const double& rho0, const double& eta, const double& mu, const double& rho, const double& eps) const;
  double compute_dphideps_rho(const double& eta, const double& eps) const;
  double compute_dP2deps_rho(const double& phi, const double& dphideps_rho, const double& rho, const double& eps) const;
  double compute_dP4deps_rho(const double& phi, const double& dphideps_rho, const double& eta, const double& rho, const double& eps) const;

  // We also want the equivalent functions for individual calculations.
  std::tuple<Scalar, Scalar, Scalar> pressureAndDerivs(const Scalar massDensity,
                                                       const Scalar specificThermalEnergy) const;

  Scalar temperature(const Scalar massDensity,
                     const Scalar specificThermalEnergy) const;

  Scalar specificThermalEnergy(const Scalar massDensity,
                               const Scalar temperature) const;

  Scalar specificHeat(const Scalar massDensity,
                      const Scalar temperature) const;

  Scalar soundSpeed(const Scalar massDensity,
                    const Scalar specificThermalEnergy) const;

  Scalar gamma(const Scalar massDensity,
               const Scalar specificThermalEnergy) const;

  Scalar bulkModulus(const Scalar massDensity,
                     const Scalar specificThermalEnergy) const;

  Scalar entropy(const Scalar massDensity,
                 const Scalar specificThermalEnergy) const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mEtaMinSolid, mEtaMaxSolid,
         ma, mb, mA, mB, malpha, mbeta, meps0, mepsLiquid, mepsVapor,
         mAtomicWeight, mCv, mdPdRhoMin;

  // Disallow default constructor
  TillotsonEquationOfState();

  using EquationOfState<Dimension>::mConstants;
};

}

#include "TillotsonEquationOfStateInline.hh"

#endif
