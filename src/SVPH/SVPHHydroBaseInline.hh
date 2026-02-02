namespace Spheral {

//------------------------------------------------------------------------------
// Choose whether we want to sum for mass density, or integrate the continuity
// equation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
MassDensityType
SVPHHydroBase<Dimension>::densityUpdate() const {
  return mDensityUpdate;
}

template<typename Dimension>
inline
void
SVPHHydroBase<Dimension>::
densityUpdate(const MassDensityType type) {
  mDensityUpdate = type;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the compatible energy evolution 
// algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SVPHHydroBase<Dimension>::compatibleEnergyEvolution() const {
  return mCompatibleEnergyEvolution;
}

template<typename Dimension>
inline
void
SVPHHydroBase<Dimension>::compatibleEnergyEvolution(const bool val) {
  mCompatibleEnergyEvolution = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the XSVPH algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SVPHHydroBase<Dimension>::XSVPH() const {
  return mXSVPH;
}

template<typename Dimension>
inline
void
SVPHHydroBase<Dimension>::XSVPH(const bool val) {
  mXSVPH = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the linear consistency corretions.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SVPHHydroBase<Dimension>::linearConsistent() const {
  return mLinearConsistent;
}

template<typename Dimension>
inline
void
SVPHHydroBase<Dimension>::linearConsistent(const bool val) {
  mLinearConsistent = val;
}

//------------------------------------------------------------------------------
// The fraction of centroidal motion to apply.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
SVPHHydroBase<Dimension>::fcentroidal() const {
  return mfcentroidal;
}

template<typename Dimension>
inline
void
SVPHHydroBase<Dimension>::fcentroidal(const typename Dimension::Scalar val) {
  VERIFY2(val >= 0.0 and val <= 1.0,
          "SVPHHydro range error : fcentroidal should be in the range [0,1].");
  mfcentroidal = val;
}

//------------------------------------------------------------------------------
// Access the optional min & max bounds for generating meshes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Vector&
SVPHHydroBase<Dimension>::
xmin() const {
  return mXmin;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
SVPHHydroBase<Dimension>::
xmax() const {
  return mXmax;
}

template<typename Dimension>
inline
void
SVPHHydroBase<Dimension>::
xmin(const typename Dimension::Vector& x) {
  mXmin = x;
}

template<typename Dimension>
inline
void
SVPHHydroBase<Dimension>::
xmax(const typename Dimension::Vector& x) {
  mXmax = x;
}

//------------------------------------------------------------------------------
// The interpolation kernel
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
SVPHHydroBase<Dimension>::
kernel() const {
  return mKernel;
}

//------------------------------------------------------------------------------
// The mesh.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const Mesh<Dimension>&
SVPHHydroBase<Dimension>::
mesh() const {
  return *mMeshPtr;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
A() const {
  return mA;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SVPHHydroBase<Dimension>::
B() const {
  return mB;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SVPHHydroBase<Dimension>::
gradB() const {
  return mGradB;
}

template<typename Dimension>
inline
const FieldList<Dimension, int>&
SVPHHydroBase<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
pressure() const {
  return mPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
soundSpeed() const {
  return mSoundSpeed;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
volume() const {
  return mVolume;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
massDensitySum() const {
  return mMassDensitySum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SVPHHydroBase<Dimension>::
XSVPHDeltaV() const {
  return mXSVPHDeltaV;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SVPHHydroBase<Dimension>::
DxDt() const {
  return mDxDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SVPHHydroBase<Dimension>::
DvDt() const {
  return mDvDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
DmassDensityDt() const {
  return mDmassDensityDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
DspecificThermalEnergyDt() const {
  return mDspecificThermalEnergyDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SVPHHydroBase<Dimension>::
DvDx() const {
  return mDvDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SVPHHydroBase<Dimension>::
internalDvDx() const {
  return mInternalDvDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, std::vector<typename Dimension::Vector> >&
SVPHHydroBase<Dimension>::
pairAccelerations() const {
  return mPairAccelerations;
}

}
