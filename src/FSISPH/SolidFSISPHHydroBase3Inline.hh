namespace Spheral {


//------------------------------------------------------------------------------
// Access the main kernel used for (A)SPH field estimates.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
SolidFSISPHHydroBase<Dimension>::
kernel() const {
  return mKernel;
}

//------------------------------------------------------------------------------
// The object defining how smoothing scales are evolved.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const SmoothingScaleBase<Dimension>&
SolidFSISPHHydroBase<Dimension>::
smoothingScaleMethod() const {
  return mSmoothingScaleMethod;
}

//------------------------------------------------------------------------------
// Ref to the slide surface obj
//------------------------------------------------------------------------------
template<typename Dimension>
inline
SlideSurface<Dimension>&
SolidFSISPHHydroBase<Dimension>::slideSurface() const {
  return mSlideSurface;
}

//------------------------------------------------------------------------------
// Choose whether we want to sum for mass density, or integrate the continuity
// equation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FSIMassDensityMethod
SolidFSISPHHydroBase<Dimension>::densityUpdate() const {
  return mDensityUpdate;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
densityUpdate(FSIMassDensityMethod type) {
  mDensityUpdate = type;
}

//------------------------------------------------------------------------------
// Choose how we want to update the H tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
HEvolutionType
SolidFSISPHHydroBase<Dimension>::HEvolution() const {
  return mHEvolution;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
HEvolution(HEvolutionType type) {
  mHEvolution = type;
}

//------------------------------------------------------------------------------
// return our interface method
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
interfaceMethod(InterfaceMethod x) {
  mInterfaceMethod = x;
}
template<typename Dimension>
inline
InterfaceMethod
SolidFSISPHHydroBase<Dimension>::
interfaceMethod() const {
  return mInterfaceMethod;
}

//------------------------------------------------------------------------------
// return our kernel method
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
kernelAveragingMethod(KernelAveragingMethod x) {
  mKernelAveragingMethod = x;
}
template<typename Dimension>
inline
KernelAveragingMethod
SolidFSISPHHydroBase<Dimension>::
kernelAveragingMethod() const {
  return mKernelAveragingMethod;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the compatible energy evolution 
// algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SolidFSISPHHydroBase<Dimension>::compatibleEnergyEvolution() const {
  return mCompatibleEnergyEvolution;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::compatibleEnergyEvolution(bool val) {
  mCompatibleEnergyEvolution = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're evolving total or specific energy
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SolidFSISPHHydroBase<Dimension>::evolveTotalEnergy() const {
  return mEvolveTotalEnergy;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::evolveTotalEnergy(bool val) {
  mEvolveTotalEnergy = val;
}

//------------------------------------------------------------------------------
// Access the flag controlling linear correct velocity gradient.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SolidFSISPHHydroBase<Dimension>::linearCorrectGradients() const {
  return mLinearCorrectGradients;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::linearCorrectGradients(bool val) {
  mLinearCorrectGradients = val;
}

//------------------------------------------------------------------------------
// plane strain (instead of plane stress) mode for 1D and 2D problems
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SolidFSISPHHydroBase<Dimension>::planeStrain() const {
  return mPlaneStrain;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::planeStrain(bool val) {
  mPlaneStrain = val;
}


//------------------------------------------------------------------------------
// switch to turn on density sum for different nodeLists
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
applySelectSumDensity(bool x) {
  mApplySelectDensitySum = x;
}
template<typename Dimension>
inline
bool
SolidFSISPHHydroBase<Dimension>::
applySelectSumDensity() const {
  return mApplySelectDensitySum;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
sumDensityNodeLists(std::vector<int> x) {
  mSumDensityNodeLists = x;
}
template<typename Dimension>
inline
std::vector<int>
SolidFSISPHHydroBase<Dimension>::
sumDensityNodeLists() const {
  return mSumDensityNodeLists;
}

//------------------------------------------------------------------------------
// set/get fsi coeffs
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
surfaceForceCoefficient(double x) {
  mSurfaceForceCoefficient = x;
}
template<typename Dimension>
inline
double
SolidFSISPHHydroBase<Dimension>::
surfaceForceCoefficient() const {
  return mSurfaceForceCoefficient;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
densityStabilizationCoefficient(double x) {
  mDensityStabilizationCoefficient = x;
}
template<typename Dimension>
inline
double
SolidFSISPHHydroBase<Dimension>::
densityStabilizationCoefficient() const {
  return mDensityStabilizationCoefficient;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
specificThermalEnergyDiffusionCoefficient(double x) {
  mSpecificThermalEnergyDiffusionCoefficient = x;
}
template<typename Dimension>
inline
double
SolidFSISPHHydroBase<Dimension>::
specificThermalEnergyDiffusionCoefficient() const {
  return mSpecificThermalEnergyDiffusionCoefficient;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
xsphCoefficient(double x) {
  mXSPHCoefficient = x;
}
template<typename Dimension>
inline
double
SolidFSISPHHydroBase<Dimension>::
xsphCoefficient() const {
  return mXSPHCoefficient;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
interfacePmin(double x) {
  mInterfacePmin = x;
}
template<typename Dimension>
inline
double
SolidFSISPHHydroBase<Dimension>::
interfacePmin() const {
  return mInterfacePmin;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
interfaceNeighborAngleThreshold(double x) {
  mInterfaceNeighborAngleThreshold = x;
}
template<typename Dimension>
inline
double
SolidFSISPHHydroBase<Dimension>::
interfaceNeighborAngleThreshold() const {
  return mInterfaceNeighborAngleThreshold;
}

//------------------------------------------------------------------------------
// Parameter to determine the magnitude of the tensile small scale correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
SolidFSISPHHydroBase<Dimension>::
epsilonTensile() const {
  return mEpsTensile;
}

template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
epsilonTensile(typename Dimension::Scalar val) {
  mEpsTensile = val;
}

//------------------------------------------------------------------------------
// Parameter to set the exponent used in the tensile small scale correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
SolidFSISPHHydroBase<Dimension>::
nTensile() const {
  return mnTensile;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
nTensile(typename Dimension::Scalar val) {
  mnTensile = val;
}

//------------------------------------------------------------------------------
// Access the optional min & max bounds for generating meshes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Vector&
SolidFSISPHHydroBase<Dimension>::
xmin() const {
  return mxmin;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
SolidFSISPHHydroBase<Dimension>::
xmax() const {
  return mxmax;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
xmin(const typename Dimension::Vector& x) {
  mxmin = x;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
xmax(const typename Dimension::Vector& x) {
  mxmax = x;
}

//------------------------------------------------------------------------------
// Get methods for FieldList refs
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
pairAccelerations() const {
  return mPairAccelerations;
}

template<typename Dimension>
inline
const std::vector<typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
pairDepsDt() const {
  return mPairDepsDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, int>&
SolidFSISPHHydroBase<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
pressure() const {
  return mPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
damagedPressure() const {
  return mDamagedPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
soundSpeed() const {
  return mSoundSpeed;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
bulkModulus() const {
  return mBulkModulus;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
shearModulus() const {
  return mShearModulus;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
yieldStrength() const {
  return mYieldStrength;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
plasticStrain0() const {
  return mPlasticStrain0;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
inverseEquivalentDeviatoricStress() const {
  return mInverseEquivalentDeviatoricStress;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
volume() const {
  return mVolume;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
DxDt() const {
  return mDxDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
XSPHDeltaV() const {
  return mXSPHDeltaV;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
XSPHWeightSum() const {
  return mXSPHWeightSum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
DvDt() const {
  return mDvDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
DmassDensityDt() const {
  return mDmassDensityDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
DspecificThermalEnergyDt() const {
  return mDspecificThermalEnergyDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
SolidFSISPHHydroBase<Dimension>::
DdeviatoricStressDt() const {
  return mDdeviatoricStressDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
SolidFSISPHHydroBase<Dimension>::
DHDt() const {
  return mDHDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
SolidFSISPHHydroBase<Dimension>::
Hideal() const {
  return mHideal;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
DPDx() const {
  return mDPDx;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
DepsDx() const {
  return mDepsDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SolidFSISPHHydroBase<Dimension>::
DvDx() const {
  return mDvDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SolidFSISPHHydroBase<Dimension>::
internalDvDx() const {
  return mInternalDvDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SolidFSISPHHydroBase<Dimension>::
M() const {
  return mM;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SolidFSISPHHydroBase<Dimension>::
localM() const {
  return mLocalM;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
maxViscousPressure() const {
  return mMaxViscousPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
normalization() const {
  return mNormalization;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::ThirdRankTensor>&
SolidFSISPHHydroBase<Dimension>::
DSDx() const {
  return mDSDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
weightedNeighborSum() const {
  return mWeightedNeighborSum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
SolidFSISPHHydroBase<Dimension>::
massSecondMoment() const {
  return mMassSecondMoment;
}

// template<typename Dimension>
// inline
// const FieldList<Dimension,  typename Dimension::Vector>&
// SolidFSISPHHydroBase<Dimension>::
// interfaceNormals() const {
//   return mInterfaceNormals;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension,  typename Dimension::Scalar>&
// SolidFSISPHHydroBase<Dimension>::
// interfaceFraction() const {
//   return mInterfaceFraction;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension,  typename Dimension::Scalar>&
// SolidFSISPHHydroBase<Dimension>::
// interfaceSmoothness() const {
//   return mInterfaceSmoothness;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension,  typename Dimension::Vector>&
// SolidFSISPHHydroBase<Dimension>::
// newInterfaceNormals() const {
//   return mNewInterfaceNormals;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension,  typename Dimension::Vector>&
// SolidFSISPHHydroBase<Dimension>::
// smoothedInterfaceNormals() const {
//   return mSmoothedInterfaceNormals;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension,  typename Dimension::Scalar>&
// SolidFSISPHHydroBase<Dimension>::
// newInterfaceFraction() const {
//   return mNewInterfaceFraction;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension,  typename Dimension::Scalar>&
// SolidFSISPHHydroBase<Dimension>::
// newInterfaceSmoothness() const {
//   return mNewInterfaceSmoothness;
// }


template<typename Dimension>
inline
const FieldList<Dimension,  int>&
SolidFSISPHHydroBase<Dimension>::
interfaceFlags() const {
  return mInterfaceFlags;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
interfaceAreaVectors() const {
  return mInterfaceAreaVectors;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
interfaceNormals() const {
  return mInterfaceNormals;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
interfaceSmoothness() const {
  return mInterfaceSmoothness;
}

template<typename Dimension>
inline
const FieldList<Dimension,  int>&
SolidFSISPHHydroBase<Dimension>::
newInterfaceFlags() const {
  return mNewInterfaceFlags;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
newInterfaceAreaVectors() const {
  return mNewInterfaceAreaVectors;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
newInterfaceNormals() const {
  return mNewInterfaceNormals;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
interfaceSmoothnessNormalization() const {
  return mInterfaceSmoothnessNormalization;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
interfaceFraction() const {
  return mInterfaceFraction;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
newInterfaceSmoothness() const {
  return mNewInterfaceSmoothness;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
interfaceAngles() const {
  return mInterfaceAngles;
}

}
