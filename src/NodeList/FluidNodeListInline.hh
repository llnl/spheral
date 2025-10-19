#include "Field/Field.hh"
#include "Utilities/SpheralFunctions.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Mass density per node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, typename Dimension::Scalar>&
FluidNodeList<Dimension>::massDensity() {
  REQUIRE(mMassDensity.nodeListPtr() == this);
  return mMassDensity;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
FluidNodeList<Dimension>::massDensity() const {
  REQUIRE(mMassDensity.nodeListPtr() == this);
  return mMassDensity;
}

//------------------------------------------------------------------------------
// Specific thermal energy per node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, typename Dimension::Scalar>&
FluidNodeList<Dimension>::specificThermalEnergy() {
  REQUIRE(mSpecificThermalEnergy.nodeListPtr() == this);
  return mSpecificThermalEnergy;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
FluidNodeList<Dimension>::specificThermalEnergy() const {
  REQUIRE(mSpecificThermalEnergy.nodeListPtr() == this);
  return mSpecificThermalEnergy;
}

//------------------------------------------------------------------------------
// Access the equation of state.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const EquationOfState<Dimension>&
FluidNodeList<Dimension>::equationOfState() const {
  return *mEosPtr;
}

template<typename Dimension>
inline
void
FluidNodeList<Dimension>::equationOfState(const EquationOfState<Dimension>& eos) {
  mEosPtr = &eos;
}

//------------------------------------------------------------------------------
// Min/max mass densities when we're time integrating.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
FluidNodeList<Dimension>::
rhoMin() const {
  return mRhoMin;
}

template<typename Dimension>
inline
void
FluidNodeList<Dimension>::
rhoMin(typename Dimension::Scalar x) {
  mRhoMin = x;
}

template<typename Dimension>
inline
typename Dimension::Scalar
FluidNodeList<Dimension>::
rhoMax() const {
  return mRhoMax;
}

template<typename Dimension>
inline
void
FluidNodeList<Dimension>::
rhoMax(typename Dimension::Scalar x) {
  mRhoMax = x;
}


//------------------------------------------------------------------------------
// allow ALE schemes in MFV
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
FluidNodeList<Dimension>::
allowALE() const {
  return mAllowALE;
}

template<typename Dimension>
inline
void
FluidNodeList<Dimension>::
allowALE(bool x) {
  mAllowALE = x;
}

}
