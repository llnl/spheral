namespace Spheral {

//------------------------------------------------------------------------------
// Access the CFL safety criteria.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
GenericRiemannHydro<Dimension>::cfl() const {
  return mCfl;
}

template<typename Dimension>
inline
void
GenericRiemannHydro<Dimension>::
cfl(typename Dimension::Scalar cfl) {
  mCfl = cfl;
}

//------------------------------------------------------------------------------
// Flag to set whether or not to use the magnitude of the velocity to set the
// timestep.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GenericRiemannHydro<Dimension>::useVelocityMagnitudeForDt() const {
  return mUseVelocityMagnitudeForDt;
}

template<typename Dimension>
inline
void
GenericRiemannHydro<Dimension>::
useVelocityMagnitudeForDt(bool x) {
  mUseVelocityMagnitudeForDt = x;
}

//------------------------------------------------------------------------------
// Return the master neighboring statistics.
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// int
// GenericRiemannHydro<Dimension>::minMasterNeighbor() const {
//   return mMinMasterNeighbor;
// }

// template<typename Dimension>
// inline
// int
// GenericRiemannHydro<Dimension>::maxMasterNeighbor() const {
//   return mMaxMasterNeighbor;
// }

// template<typename Dimension>
// inline
// double
// GenericRiemannHydro<Dimension>::averageMasterNeighbor() const {
//   return double(mSumMasterNeighbor)/(mNormMasterNeighbor + FLT_MIN);
// }

//------------------------------------------------------------------------------
// Return the coarse neighboring statistics.
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// int
// GenericRiemannHydro<Dimension>::minCoarseNeighbor() const {
//   return mMinCoarseNeighbor;
// }

// template<typename Dimension>
// inline
// int
// GenericRiemannHydro<Dimension>::maxCoarseNeighbor() const {
//   return mMaxCoarseNeighbor;
// }

// template<typename Dimension>
// inline
// double
// GenericRiemannHydro<Dimension>::averageCoarseNeighbor() const {
//   return double(mSumCoarseNeighbor)/(mNormCoarseNeighbor + FLT_MIN);
// }

//------------------------------------------------------------------------------
// Return the refine neighboring statistics.
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// int
// GenericRiemannHydro<Dimension>::minRefineNeighbor() const {
//   return mMinRefineNeighbor;
// }

// template<typename Dimension>
// inline
// int
// GenericRiemannHydro<Dimension>::maxRefineNeighbor() const {
//   return mMaxRefineNeighbor;
// }

// template<typename Dimension>
// inline
// double
// GenericRiemannHydro<Dimension>::averageRefineNeighbor() const {
//   return double(mSumRefineNeighbor)/(mNormRefineNeighbor + FLT_MIN);
// }

//------------------------------------------------------------------------------
// Return the actual neighboring statistics.
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// int
// GenericRiemannHydro<Dimension>::minActualNeighbor() const {
//   return mMinActualNeighbor;
// }

// template<typename Dimension>
// inline
// int
// GenericRiemannHydro<Dimension>::maxActualNeighbor() const {
//   return mMaxActualNeighbor;
// }

// template<typename Dimension>
// inline
// double
// GenericRiemannHydro<Dimension>::averageActualNeighbor() const {
//   return double(mSumActualNeighbor)/(mNormActualNeighbor + FLT_MIN);
// }

//------------------------------------------------------------------------------
// Update the master neighboring statistics.
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// void
// GenericRiemannHydro<Dimension>::
// updateMasterNeighborStats(int numMaster) const {
//   if (numMaster > 0) {
//     mMinMasterNeighbor = std::min(mMinMasterNeighbor, numMaster);
//     mMaxMasterNeighbor = std::max(mMaxMasterNeighbor, numMaster);
//     mSumMasterNeighbor += numMaster;
//     mNormMasterNeighbor++;
//   }
// }

//------------------------------------------------------------------------------
// Update the coarse neighboring statistics.
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// void
// GenericRiemannHydro<Dimension>::
// updateCoarseNeighborStats(int numNeighbor) const {
//   if (numNeighbor > 0) {
//     mMinCoarseNeighbor = std::min(mMinCoarseNeighbor, numNeighbor);
//     mMaxCoarseNeighbor = std::max(mMaxCoarseNeighbor, numNeighbor);
//     mSumCoarseNeighbor += numNeighbor;
//     mNormCoarseNeighbor++;
//   }
// }

//------------------------------------------------------------------------------
// Update the refine neighboring statistics.
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// void
// GenericRiemannHydro<Dimension>::
// updateRefineNeighborStats(int numNeighbor) const {
//   if (numNeighbor > 0) {
//     mMinRefineNeighbor = std::min(mMinRefineNeighbor, numNeighbor);
//     mMaxRefineNeighbor = std::max(mMaxRefineNeighbor, numNeighbor);
//     mSumRefineNeighbor += numNeighbor;
//     mNormRefineNeighbor++;
//   }
// }

//------------------------------------------------------------------------------
// Update the actual neighboring statistics.
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// void
// GenericRiemannHydro<Dimension>::
// updateActualNeighborStats(int numNeighbor) const {
//   if (numNeighbor > 0) {
//     mMinActualNeighbor = std::min(mMinActualNeighbor, numNeighbor);
//     mMaxActualNeighbor = std::max(mMaxActualNeighbor, numNeighbor);
//     mSumActualNeighbor += numNeighbor;
//     mNormActualNeighbor++;
//   }
// }

//------------------------------------------------------------------------------
// get out reimann solver obj
//------------------------------------------------------------------------------
template<typename Dimension>
inline
RiemannSolverBase<Dimension>&
GenericRiemannHydro<Dimension>::riemannSolver() const {
  return mRiemannSolver;
}


//------------------------------------------------------------------------------
// set/get gradient type
//------------------------------------------------------------------------------
template<typename Dimension>
inline
GradientType
GenericRiemannHydro<Dimension>::
gradientType() const {
  return mGradientType;
}

template<typename Dimension>
inline
void
GenericRiemannHydro<Dimension>::
gradientType(GradientType x) {
  mGradientType=x;
}

//------------------------------------------------------------------------------
// set/get for how we evolve our density
//------------------------------------------------------------------------------
template<typename Dimension>
inline
MassDensityType
GenericRiemannHydro<Dimension>::densityUpdate() const {
  return mDensityUpdate;
}

template<typename Dimension>
inline
void
GenericRiemannHydro<Dimension>::
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
GenericRiemannHydro<Dimension>::compatibleEnergyEvolution() const {
  return mCompatibleEnergyEvolution;
}

template<typename Dimension>
inline
void
GenericRiemannHydro<Dimension>::compatibleEnergyEvolution(bool val) {
  mCompatibleEnergyEvolution = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're evolving total or specific energy
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GenericRiemannHydro<Dimension>::evolveTotalEnergy() const {
  return mEvolveTotalEnergy;
}

template<typename Dimension>
inline
void
GenericRiemannHydro<Dimension>::evolveTotalEnergy(bool val) {
  mEvolveTotalEnergy = val;
}


//------------------------------------------------------------------------------
// Access the flag determining if we're using the XSPH algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GenericRiemannHydro<Dimension>::XSPH() const {
  return mXSPH;
}

template<typename Dimension>
inline
void
GenericRiemannHydro<Dimension>::XSPH(bool val) {
  mXSPH = val;
}

//------------------------------------------------------------------------------
// Access the flag controlling linear correct velocity gradient.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GenericRiemannHydro<Dimension>::correctVelocityGradient() const {
  return mCorrectVelocityGradient;
}

template<typename Dimension>
inline
void
GenericRiemannHydro<Dimension>::correctVelocityGradient(bool val) {
  mCorrectVelocityGradient = val;
}


//------------------------------------------------------------------------------
// Parameter to determine the magnitude of the tensile small scale correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
GenericRiemannHydro<Dimension>::
epsilonTensile() const {
  return mEpsTensile;
}

template<typename Dimension>
void
GenericRiemannHydro<Dimension>::
epsilonTensile(typename Dimension::Scalar val) {
  mEpsTensile = val;
}

//------------------------------------------------------------------------------
// Parameter to set the exponent used in the tensile small scale correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
GenericRiemannHydro<Dimension>::
nTensile() const {
  return mnTensile;
}

template<typename Dimension>
inline
void
GenericRiemannHydro<Dimension>::
nTensile(typename Dimension::Scalar val) {
  mnTensile = val;
}

//------------------------------------------------------------------------------
// Access the optional min & max bounds for generating meshes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Vector&
GenericRiemannHydro<Dimension>::
xmin() const {
  return mxmin;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
GenericRiemannHydro<Dimension>::
xmax() const {
  return mxmax;
}

template<typename Dimension>
inline
void
GenericRiemannHydro<Dimension>::
xmin(const typename Dimension::Vector& x) {
  mxmin = x;
}

template<typename Dimension>
inline
void
GenericRiemannHydro<Dimension>::
xmax(const typename Dimension::Vector& x) {
  mxmax = x;
}

//------------------------------------------------------------------------------
// Access the main kernel used for (A)SPH field estimates.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
GenericRiemannHydro<Dimension>::
kernel() const {
  return mKernel;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, int>&
GenericRiemannHydro<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
GenericRiemannHydro<Dimension>::
volume() const {
  return mVolume;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
GenericRiemannHydro<Dimension>::
pressure() const {
  return mPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
GenericRiemannHydro<Dimension>::
soundSpeed() const {
  return mSoundSpeed;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
GenericRiemannHydro<Dimension>::
normalization() const {
  return mNormalization;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
GenericRiemannHydro<Dimension>::
XSPHWeightSum() const {
  return mXSPHWeightSum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
GenericRiemannHydro<Dimension>::
XSPHDeltaV() const {
  return mXSPHDeltaV;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
GenericRiemannHydro<Dimension>::
M() const {
  return mM;
}

// template<typename Dimension>
// inline
// const FieldList<Dimension, typename Dimension::Tensor>&
// GenericRiemannHydro<Dimension>::
// localM() const {
//   return mLocalM;
// }

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
GenericRiemannHydro<Dimension>::
DxDt() const {
  return mDxDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
GenericRiemannHydro<Dimension>::
DvDt() const {
  return mDvDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
GenericRiemannHydro<Dimension>::
DspecificThermalEnergyDt() const {
  return mDspecificThermalEnergyDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
GenericRiemannHydro<Dimension>::
DvDx() const {
  return mDvDx;
}

// template<typename Dimension>
// inline
// const FieldList<Dimension, typename Dimension::Tensor>&
// GenericRiemannHydro<Dimension>::
// internalDvDx() const {
//   return mInternalDvDx;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension, typename Dimension::Vector>&
// GenericRiemannHydro<Dimension>::
// DpDx() const {
//   return mDpDx;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension, typename Dimension::Vector>&
// GenericRiemannHydro<Dimension>::
// DpDxRaw() const {
//   return mDpDxRaw;
// }

// template<typename Dimension>
// inline
// const FieldList<Dimension, typename Dimension::Tensor>&
// GenericRiemannHydro<Dimension>::
// DvDxRaw() const {
//   return mDvDxRaw;
// }



// template<typename Dimension>
// inline
// const FieldList<Dimension, typename Dimension::Vector>&
// GenericRiemannHydro<Dimension>::
// DrhoDx() const {
//   return mDrhoDx;
// }

template<typename Dimension>
inline
const typename GenericRiemannHydro<Dimension>::PairAccelerationsType&
GenericRiemannHydro<Dimension>::
pairAccelerations() const {
  VERIFY2(mPairAccelerationsPtr, "GenericRiemannHydro ERROR: attempt to access uninitialized pairAccelerations");
  return *mPairAccelerationsPtr;
}

template<typename Dimension>
inline
const typename GenericRiemannHydro<Dimension>::PairWorkType&
GenericRiemannHydro<Dimension>::
pairDepsDt() const {
  VERIFY2(mPairDepsDtPtr, "GenericRiemannHydro ERROR: attempt to access uninitialized pairDepsDt");
  return *mPairDepsDtPtr;
}

template<typename Dimension>
inline
void
GenericRiemannHydro<Dimension>::
specificThermalEnergyDiffusionCoefficient(const typename Dimension::Scalar x) {
  mSpecificThermalEnergyDiffusionCoefficient = x;
}
template<typename Dimension>
inline
typename Dimension::Scalar
GenericRiemannHydro<Dimension>::
specificThermalEnergyDiffusionCoefficient() const {
  return mSpecificThermalEnergyDiffusionCoefficient;
}


template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
GenericRiemannHydro<Dimension>::
DrhoDx() const {
  return mDrhoDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
GenericRiemannHydro<Dimension>::
riemannDpDx() const {
  return mRiemannDpDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
GenericRiemannHydro<Dimension>::
riemannDvDx() const {
  return mRiemannDvDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
GenericRiemannHydro<Dimension>::
newRiemannDpDx() const {
  return mNewRiemannDpDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
GenericRiemannHydro<Dimension>::
newRiemannDvDx() const {
  return mNewRiemannDvDx;
}

}
