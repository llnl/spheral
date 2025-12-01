//---------------------------------Spheral++----------------------------------//
// SteinbergGuinanLundStrength -- Implements the Steinberg, Guinan, & Lund
// rate dependent strength model.
//
//   Steinberg, D.J. & Lund, C.M. (1988), J. App. Phys., 65, 1528.
//
// Created by JMO, Wed Sep 8 15:18:44 2004
//----------------------------------------------------------------------------//
#include "SteinbergGuinanLundStrength.hh"
#include "Utilities/newtonRaphson.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Functor for use looking up the thermal yeild via our Newton-Raphson function.
//------------------------------------------------------------------------------
class ThermalYieldFunctor {
 public:
  ThermalYieldFunctor(const double epspdot,
                      const double C1,
                      const double C2,
                      const double UK,
                      const double T,
                      const double YP):
    mepspdot(epspdot),
    mC1(C1),
    mC2(C2),
    mUK(UK),
    mT(T),
    mYP(YP) {
      REQUIRE(mepspdot >= 0.0);
      REQUIRE(distinctlyGreaterThan(mC1, 0.0));
      REQUIRE(distinctlyGreaterThan(mT, 0.0));
      REQUIRE(distinctlyGreaterThan(mYP, 0.0));
    }
  ~ThermalYieldFunctor() {}

  std::pair<double, double> operator()(const double x) const {
    REQUIRE(x > 0.0);
    if (fuzzyEqual(mepspdot, 0.0)) {
      return std::pair<double, double>(0.0, 0.0);
    } else {
      REQUIRE(distinctlyGreaterThan(1.0/mepspdot, mC2/x));
      const double A = 2.0*mUK/mT;
      const double B = 1.0 - x/mYP;
      const double C = 1.0/mepspdot - mC2/x;
      CHECK(distinctlyGreaterThan(C, 0.0));
      return std::pair<double, double>(A*B*B - log(C) - log(mC1),
                                       -2.0*A/mYP*B + mC2/C);
    }
  }

 private:
  double mepspdot, mC1, mC2, mUK, mT, mYP;
};

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SteinbergGuinanLundStrength<Dimension>::
SteinbergGuinanLundStrength(const SolidEquationOfState<Dimension>& eos,
                            const double G0,     
                            const double A,      
                            const double B,      
                            const double Y0,     
                            const double Ymax,
                            const double Yp,
                            const double beta,
                            const double gamma0, 
                            const double nhard,
                            const double C1,
                            const double C2,
                            const double UK,
                            const double YP,
                            const double YTmax,
                            const NinthOrderPolynomialFit& coldEnergyFit,
                            const NinthOrderPolynomialFit& meltEnergyFit):
  SteinbergGuinanStrength<Dimension>(eos,
                                                G0,
                                                A,
                                                B,
                                                Y0,
                                                Ymax,
                                                Yp,
                                                beta,
                                                gamma0,
                                                nhard,
                                                coldEnergyFit,
                                                meltEnergyFit),
  mC1(C1),
  mC2(C2),
  mUK(UK),
  mYP(YP),
  mYTmax(YTmax) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SteinbergGuinanLundStrength<Dimension>::
~SteinbergGuinanLundStrength() {
}

//------------------------------------------------------------------------------
// Compute the yield strength.
//------------------------------------------------------------------------------
template<typename Dimension>
double
SteinbergGuinanLundStrength<Dimension>::
yieldStrength(const double density,
              const double specificThermalEnergy,
              const double pressure,
              const double plasticStrain,
              const double plasticStrainRate) const {

  // Build a functor to find the thermal yeild.
  const double T = this->computeTemperature(density, specificThermalEnergy);
  ThermalYieldFunctor func(plasticStrainRate,
                           mC1,
                           mC2,
                           mUK,
                           T,
                           mYP);

  // Determine a range of yield space that encompasses the root we seek.
  const double minx = max(1.0001*mC2*plasticStrainRate, 1.0e-10);
  CHECK(distinctlyGreaterThan(minx, mC2*plasticStrainRate));
  double x1 = max(minx, 1e5);
  double x2 = x1;
  int i = 0;
  while (i < 30 && (func(x1).first * func(x2).first > 0.0)) {
    ++i;
    x1 = max(minx, 0.1*x1);
    x2 *= 10.0;
  }
  CHECK2(func(x1).first * func(x2).first > 0.0,
         "SteinbergGuinanLundStrength::yieldStrength: unable to find range for yield lookup : \n"
         << density << " "
         << specificThermalEnergy << " "
         << T << " "
         << pressure << " "
         << plasticStrain << " "
         << plasticStrainRate << " "
         << minx << " "
         << x1 << " "
         << x2 << " "
         << func(x1).first << " "
         << func(x2).first)

  // Look up the thermal yield.
  const double YT = newtonRaphson(func, x1, x2);

  // Get the athermal yield.
  const double yieldAthermal = SteinbergGuinanStrength<Dimension>::yieldStrength(density,
                                                                                 specificThermalEnergy,
                                                                                 pressure,
                                                                                 plasticStrain,
                                                                                 plasticStrainRate);

  // Get the shear modulus.
  const double G = this->shearModulus(density, specificThermalEnergy, pressure);

  // Now combine these yields for the final answer.
  CHECK(distinctlyGreaterThan(this->G0(), 0.0));
  return min(YT, mYTmax)*G/this->G0() + yieldAthermal;
}

//------------------------------------------------------------------------------
// Access the strength parameters.
//------------------------------------------------------------------------------
template<typename Dimension>
double
SteinbergGuinanLundStrength<Dimension>::
C1() const {
  return mC1;
}

template<typename Dimension>
double
SteinbergGuinanLundStrength<Dimension>::
C2() const {
  return mC2;
}

template<typename Dimension>
double
SteinbergGuinanLundStrength<Dimension>::
UK() const {
  return mUK;
}

template<typename Dimension>
double
SteinbergGuinanLundStrength<Dimension>::
YP() const {
  return mYP;
}

template<typename Dimension>
double
SteinbergGuinanLundStrength<Dimension>::
YTmax() const {
  return mYTmax;
}

}
