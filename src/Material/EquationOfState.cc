//---------------------------------Spheral++----------------------------------//
// EquationOfState -- Abstract base class for the equation of state classes.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//

#include "EquationOfState.hh"
#include "Geometry/Dimension.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/bisectRoot.hh"

namespace Spheral {

namespace {

template<typename Dimension>
struct Pfunctor {
  const EquationOfState<Dimension>& mEOS;
  NodeList<Dimension> mNodes;
  mutable Field<Dimension, double> mRho0, mEps, mP;
  double mP0;

  Pfunctor(const EquationOfState<Dimension>& eos, const double rho0, const double P0):
    mEOS(eos),
    mNodes("__dummynodes__", 1, 0),
    mRho0("rho", mNodes, rho0),
    mEps("eps", mNodes),
    mP("P", mNodes),
    mP0(P0) {}

  double operator()(const double eps) const {
    mEps[0] = eps;
    mEOS.setPressure(mP, mRho0, mEps);
    return mP[0] - mP0;
  }
};

}

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
EquationOfState<Dimension>::EquationOfState(const PhysicalConstants& constants,
                                            const double minimumPressure,
                                            const double maximumPressure,
                                            const MaterialPressureMinType minPressureType,
                                            const double externalPressure):
  mConstants(constants),
  mMinimumPressure(minimumPressure),
  mMaximumPressure(maximumPressure),
  mExternalPressure(externalPressure),
  mMinPressureType(minPressureType) {
  REQUIRE(mMinimumPressure <= mMaximumPressure);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
EquationOfState<Dimension>::~EquationOfState() {
}

//------------------------------------------------------------------------------
// Look up an energy that gives the requested pressure at the specified density.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
EquationOfState<Dimension>::
specificThermalEnergyForPressure(const typename Dimension::Scalar Ptarget,
                                 const typename Dimension::Scalar rho,
                                 const typename Dimension::Scalar epsMin,
                                 const typename Dimension::Scalar epsMax,
                                 const typename Dimension::Scalar epsTol,
                                 const typename Dimension::Scalar Ptol,
                                 const unsigned maxIterations,
                                 const bool verbose) const {
  const Pfunctor<Dimension> pfunc(*this, rho, Ptarget);
  return bisectRoot(pfunc, epsMin, epsMax, epsTol, Ptol, maxIterations, verbose);
}

}
