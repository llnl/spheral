namespace Spheral {

//------------------------------------------------------------------------------
// Choose whether we want to sum for mass density, or integrate the continuity
// equation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
MassDensityType
SVPHFacetedHydroBase<Dimension>::densityUpdate() const {
  return mDensityUpdate;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::
densityUpdate(MassDensityType type) {
  mDensityUpdate = type;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the compatible energy evolution 
// algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SVPHFacetedHydroBase<Dimension>::compatibleEnergyEvolution() const {
  return mCompatibleEnergyEvolution;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::compatibleEnergyEvolution(bool val) {
  mCompatibleEnergyEvolution = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the XSVPH algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SVPHFacetedHydroBase<Dimension>::XSVPH() const {
  return mXSVPH;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::XSVPH(bool val) {
  mXSVPH = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the linear consistency corretions.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SVPHFacetedHydroBase<Dimension>::linearConsistent() const {
  return mLinearConsistent;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::linearConsistent(bool val) {
  mLinearConsistent = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're generating void points.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SVPHFacetedHydroBase<Dimension>::generateVoid() const {
  return mGenerateVoid;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::generateVoid(bool val) {
  mGenerateVoid = val;
}

//------------------------------------------------------------------------------
// The fraction of centroidal motion to apply.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
SVPHFacetedHydroBase<Dimension>::fcentroidal() const {
  return mfcentroidal;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::fcentroidal(typename Dimension::Scalar val) {
  VERIFY2(val >= 0.0 and val <= 1.0,
          "SVPHFacetedHydro range error : fcentroidal should be in the range [0,1].");
  mfcentroidal = val;
}

//------------------------------------------------------------------------------
// The fraction of local cell pressure to use.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
SVPHFacetedHydroBase<Dimension>::fcellPressure() const {
  return mfcellPressure;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::fcellPressure(typename Dimension::Scalar val) {
  VERIFY2(val >= 0.0 and val <= 1.0,
          "SVPHFacetedHydro range error : fcellPressure should be in the range [0,1].");
  mfcellPressure = val;
}

//------------------------------------------------------------------------------
// Access the optional min & max bounds for generating meshes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Vector&
SVPHFacetedHydroBase<Dimension>::
xmin() const {
  return mXmin;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
SVPHFacetedHydroBase<Dimension>::
xmax() const {
  return mXmax;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::
xmin(const typename Dimension::Vector& x) {
  mXmin = x;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::
xmax(const typename Dimension::Vector& x) {
  mXmax = x;
}

//------------------------------------------------------------------------------
// Access the main kernel used for (A)SPH field estimates.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
SVPHFacetedHydroBase<Dimension>::
kernel() const {
  return mKernel;
}

//------------------------------------------------------------------------------
// The mesh.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const Mesh<Dimension>&
SVPHFacetedHydroBase<Dimension>::
mesh() const {
  return *mMeshPtr;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// const FieldList<Dimension, typename Dimension::Scalar>&
// SVPHFacetedHydroBase<Dimension>::
// A() const {
//   return mA;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension, typename Dimension::Vector>&
// SVPHFacetedHydroBase<Dimension>::
// B() const {
//   return mB;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension, typename Dimension::Tensor>&
// SVPHFacetedHydroBase<Dimension>::
// gradB() const {
//   return mGradB;
// }

template<typename Dimension>
inline
const FieldList<Dimension, int>&
SVPHFacetedHydroBase<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
pressure() const {
  return mPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
cellPressure() const {
  return mCellPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
soundSpeed() const {
  return mSoundSpeed;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
volume() const {
  return mVolume;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
massDensitySum() const {
  return mMassDensitySum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SVPHFacetedHydroBase<Dimension>::
XSVPHDeltaV() const {
  return mXSVPHDeltaV;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SVPHFacetedHydroBase<Dimension>::
DxDt() const {
  return mDxDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SVPHFacetedHydroBase<Dimension>::
DvDt() const {
  return mDvDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
DmassDensityDt() const {
  return mDmassDensityDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
DspecificThermalEnergyDt() const {
  return mDspecificThermalEnergyDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SVPHFacetedHydroBase<Dimension>::
DvDx() const {
  return mDvDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SVPHFacetedHydroBase<Dimension>::
internalDvDx() const {
  return mInternalDvDx;
}

// template<typename Dimension>
// inline
// const FieldList<Dimension, std::vector<typename Dimension::Scalar> >&
// SVPHFacetedHydroBase<Dimension>::
// faceMass() const {
//   return mFaceMass;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension, std::vector<typename Dimension::Vector> >&
// SVPHFacetedHydroBase<Dimension>::
// faceVelocity() const {
//   return mFaceVelocity;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension, std::vector<typename Dimension::Vector> >&
// SVPHFacetedHydroBase<Dimension>::
// faceAcceleration() const {
//   return mFaceAcceleration;
// }

template<typename Dimension>
inline
const FieldList<Dimension, std::vector<typename Dimension::Vector> >&
SVPHFacetedHydroBase<Dimension>::
faceForce() const {
  return mFaceForce;
}

}
