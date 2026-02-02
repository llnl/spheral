namespace Spheral {


//------------------------------------------------------------------------------
// set/get for mesh motion coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
MFV<Dimension>::nodeMotionCoefficient() const {
  return mNodeMotionCoefficient;
}

template<typename Dimension>
inline
void
MFV<Dimension>::
nodeMotionCoefficient(typename Dimension::Scalar x) {
  mNodeMotionCoefficient = x;
}

//------------------------------------------------------------------------------
// set/get mesh motion type
//------------------------------------------------------------------------------
template<typename Dimension>
inline
NodeMotionType
MFV<Dimension>::
nodeMotionType() const {
  return mNodeMotionType;
}

template<typename Dimension>
inline
void
MFV<Dimension>::
nodeMotionType(NodeMotionType x) {
  mNodeMotionType=x;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
MFV<Dimension>::
nodalVelocity() const {
  return mNodalVelocity;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
MFV<Dimension>::
DmassDt() const {
  return mDmassDt;
}
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
MFV<Dimension>::
DthermalEnergyDt() const {
  return mDthermalEnergyDt;
}
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
MFV<Dimension>::
DmomentumDt() const {
  return mDmomentumDt;
}
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
MFV<Dimension>::
DvolumeDt() const {
  return mDvolumeDt;
}
template<typename Dimension>
inline
const typename MFV<Dimension>::PairMassFluxType&
MFV<Dimension>::
pairMassFlux() const {
  VERIFY2(mPairMassFluxPtr, "MFV ERROR: attempt to access uninitialized pairMassFlux");
  return *mPairMassFluxPtr;
}

// template<typename Dimension>
// inline
// const FieldList<Dimension, typename Dimension::SymTensor>&
// MFV<Dimension>::
// HStretchTensor() const {
//   return mHStretchTensor;
// }
}
