namespace Spheral {


//------------------------------------------------------------------------------
// Access the main kernel used for (A)SPH field estimates.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
SolidFSISPH<Dimension>::
kernel() const {
  return mKernel;
}

//------------------------------------------------------------------------------
// Ref to the slide surface obj
//------------------------------------------------------------------------------
template<typename Dimension>
inline
SlideSurface<Dimension>&
SolidFSISPH<Dimension>::slideSurface() const {
  return mSlideSurface;
}

//------------------------------------------------------------------------------
// Choose whether we want to sum for mass density, or integrate the continuity
// equation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FSIMassDensityMethod
SolidFSISPH<Dimension>::densityUpdate() const {
  return mDensityUpdate;
}

template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::
densityUpdate(FSIMassDensityMethod type) {
  mDensityUpdate = type;
}

//------------------------------------------------------------------------------
// return our interface method
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::
interfaceMethod(InterfaceMethod x) {
  mInterfaceMethod = x;
}
template<typename Dimension>
inline
InterfaceMethod
SolidFSISPH<Dimension>::
interfaceMethod() const {
  return mInterfaceMethod;
}

//------------------------------------------------------------------------------
// return our kernel method
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::
kernelAveragingMethod(KernelAveragingMethod x) {
  mKernelAveragingMethod = x;
}
template<typename Dimension>
inline
KernelAveragingMethod
SolidFSISPH<Dimension>::
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
SolidFSISPH<Dimension>::compatibleEnergyEvolution() const {
  return mCompatibleEnergyEvolution;
}

template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::compatibleEnergyEvolution(bool val) {
  mCompatibleEnergyEvolution = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're evolving total or specific energy
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SolidFSISPH<Dimension>::evolveTotalEnergy() const {
  return mEvolveTotalEnergy;
}

template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::evolveTotalEnergy(bool val) {
  mEvolveTotalEnergy = val;
}

//------------------------------------------------------------------------------
// Access the flag controlling linear correct velocity gradient.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SolidFSISPH<Dimension>::linearCorrectGradients() const {
  return mLinearCorrectGradients;
}

template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::linearCorrectGradients(bool val) {
  mLinearCorrectGradients = val;
}

//------------------------------------------------------------------------------
// deactivate our decoupling algo for same material damaged nodes
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::
decoupleDamagedMaterial(bool x) {
  mDecoupleDamagedMaterial = x;
}
template<typename Dimension>
inline
bool
SolidFSISPH<Dimension>::
decoupleDamagedMaterial() const {
  return mDecoupleDamagedMaterial;
}

//------------------------------------------------------------------------------
// switch to turn on density sum for different nodeLists
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::
applySelectSumDensity(bool x) {
  mApplySelectDensitySum = x;
}
template<typename Dimension>
inline
bool
SolidFSISPH<Dimension>::
applySelectSumDensity() const {
  return mApplySelectDensitySum;
}

template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::
sumDensityNodeLists(std::vector<int> x) {
  mSumDensityNodeLists = x;
}
template<typename Dimension>
inline
std::vector<int>
SolidFSISPH<Dimension>::
sumDensityNodeLists() const {
  return mSumDensityNodeLists;
}

//------------------------------------------------------------------------------
// set/get fsi coeffs
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::
surfaceForceCoefficient(double x) {
  mSurfaceForceCoefficient = x;
}
template<typename Dimension>
inline
double
SolidFSISPH<Dimension>::
surfaceForceCoefficient() const {
  return mSurfaceForceCoefficient;
}

template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::
densityStabilizationCoefficient(double x) {
  mDensityStabilizationCoefficient = x;
}
template<typename Dimension>
inline
double
SolidFSISPH<Dimension>::
densityStabilizationCoefficient() const {
  return mDensityStabilizationCoefficient;
}

template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::
specificThermalEnergyDiffusionCoefficient(double x) {
  mSpecificThermalEnergyDiffusionCoefficient = x;
}
template<typename Dimension>
inline
double
SolidFSISPH<Dimension>::
specificThermalEnergyDiffusionCoefficient() const {
  return mSpecificThermalEnergyDiffusionCoefficient;
}

template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::
xsphCoefficient(double x) {
  mXSPHCoefficient = x;
}
template<typename Dimension>
inline
double
SolidFSISPH<Dimension>::
xsphCoefficient() const {
  return mXSPHCoefficient;
}

template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::
interfacePmin(double x) {
  mInterfacePmin = x;
}
template<typename Dimension>
inline
double
SolidFSISPH<Dimension>::
interfacePmin() const {
  return mInterfacePmin;
}

template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::
interfaceNeighborAngleThreshold(double x) {
  mInterfaceNeighborAngleThreshold = x;
}
template<typename Dimension>
inline
double
SolidFSISPH<Dimension>::
interfaceNeighborAngleThreshold() const {
  return mInterfaceNeighborAngleThreshold;
}

//------------------------------------------------------------------------------
// Parameter to determine the magnitude of the tensile small scale correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
SolidFSISPH<Dimension>::
epsilonTensile() const {
  return mEpsTensile;
}

template<typename Dimension>
void
SolidFSISPH<Dimension>::
epsilonTensile(typename Dimension::Scalar val) {
  mEpsTensile = val;
}

//------------------------------------------------------------------------------
// Parameter to set the exponent used in the tensile small scale correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
SolidFSISPH<Dimension>::
nTensile() const {
  return mnTensile;
}

template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::
nTensile(typename Dimension::Scalar val) {
  mnTensile = val;
}

//------------------------------------------------------------------------------
// Access the optional min & max bounds for generating meshes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Vector&
SolidFSISPH<Dimension>::
xmin() const {
  return mxmin;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
SolidFSISPH<Dimension>::
xmax() const {
  return mxmax;
}

template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::
xmin(const typename Dimension::Vector& x) {
  mxmin = x;
}

template<typename Dimension>
inline
void
SolidFSISPH<Dimension>::
xmax(const typename Dimension::Vector& x) {
  mxmax = x;
}

//------------------------------------------------------------------------------
// Get methods for FieldList refs
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename SolidFSISPH<Dimension>::PairAccelerationsType&
SolidFSISPH<Dimension>::
pairAccelerations() const {
  VERIFY2(mPairAccelerationsPtr, "SolidFSISPH ERROR: attempt to access uninitialized pairAccelerations");
  return *mPairAccelerationsPtr;
}

template<typename Dimension>
inline
const typename SolidFSISPH<Dimension>::PairWorkType&
SolidFSISPH<Dimension>::
pairDepsDt() const {
  VERIFY2(mPairDepsDtPtr, "SolidFSISPH ERROR: attempt to access uninitialized pairDepsDt");
  return *mPairDepsDtPtr;
}

template<typename Dimension>
inline
const FieldList<Dimension, int>&
SolidFSISPH<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
pressure() const {
  return mPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
damagedPressure() const {
  return mDamagedPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
soundSpeed() const {
  return mSoundSpeed;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
bulkModulus() const {
  return mBulkModulus;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
shearModulus() const {
  return mShearModulus;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
yieldStrength() const {
  return mYieldStrength;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
plasticStrain0() const {
  return mPlasticStrain0;
}

// template<typename Dimension>
// inline
// const FieldList<Dimension, typename Dimension::Scalar>&
// SolidFSISPH<Dimension>::
// inverseEquivalentDeviatoricStress() const {
//   return mInverseEquivalentDeviatoricStress;
// }

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
volume() const {
  return mVolume;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SolidFSISPH<Dimension>::
DxDt() const {
  return mDxDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SolidFSISPH<Dimension>::
XSPHDeltaV() const {
  return mXSPHDeltaV;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
XSPHWeightSum() const {
  return mXSPHWeightSum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SolidFSISPH<Dimension>::
DvDt() const {
  return mDvDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
DmassDensityDt() const {
  return mDmassDensityDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
DspecificThermalEnergyDt() const {
  return mDspecificThermalEnergyDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
SolidFSISPH<Dimension>::
DdeviatoricStressDt() const {
  return mDdeviatoricStressDt;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPH<Dimension>::
DPDx() const {
  return mDPDx;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPH<Dimension>::
DepsDx() const {
  return mDepsDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SolidFSISPH<Dimension>::
DvDx() const {
  return mDvDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SolidFSISPH<Dimension>::
internalDvDx() const {
  return mInternalDvDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SolidFSISPH<Dimension>::
M() const {
  return mM;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SolidFSISPH<Dimension>::
localM() const {
  return mLocalM;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
normalization() const {
  return mNormalization;
}

// template<typename Dimension>
// inline
// const FieldList<Dimension,  typename Dimension::Vector>&
// SolidFSISPH<Dimension>::
// interfaceNormals() const {
//   return mInterfaceNormals;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension,  typename Dimension::Scalar>&
// SolidFSISPH<Dimension>::
// interfaceFraction() const {
//   return mInterfaceFraction;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension,  typename Dimension::Scalar>&
// SolidFSISPH<Dimension>::
// interfaceSmoothness() const {
//   return mInterfaceSmoothness;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension,  typename Dimension::Vector>&
// SolidFSISPH<Dimension>::
// newInterfaceNormals() const {
//   return mNewInterfaceNormals;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension,  typename Dimension::Vector>&
// SolidFSISPH<Dimension>::
// smoothedInterfaceNormals() const {
//   return mSmoothedInterfaceNormals;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension,  typename Dimension::Scalar>&
// SolidFSISPH<Dimension>::
// newInterfaceFraction() const {
//   return mNewInterfaceFraction;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension,  typename Dimension::Scalar>&
// SolidFSISPH<Dimension>::
// newInterfaceSmoothness() const {
//   return mNewInterfaceSmoothness;
// }


template<typename Dimension>
inline
const FieldList<Dimension,  int>&
SolidFSISPH<Dimension>::
interfaceFlags() const {
  return mInterfaceFlags;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPH<Dimension>::
interfaceAreaVectors() const {
  return mInterfaceAreaVectors;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPH<Dimension>::
interfaceNormals() const {
  return mInterfaceNormals;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
interfaceSmoothness() const {
  return mInterfaceSmoothness;
}

template<typename Dimension>
inline
const FieldList<Dimension,  int>&
SolidFSISPH<Dimension>::
newInterfaceFlags() const {
  return mNewInterfaceFlags;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPH<Dimension>::
newInterfaceAreaVectors() const {
  return mNewInterfaceAreaVectors;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPH<Dimension>::
newInterfaceNormals() const {
  return mNewInterfaceNormals;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
interfaceSmoothnessNormalization() const {
  return mInterfaceSmoothnessNormalization;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
interfaceFraction() const {
  return mInterfaceFraction;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
newInterfaceSmoothness() const {
  return mNewInterfaceSmoothness;
}

template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPH<Dimension>::
interfaceAngles() const {
  return mInterfaceAngles;
}

}
