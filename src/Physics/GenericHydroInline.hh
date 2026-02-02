#include "Utilities/DBC.hh"

#include <limits.h>
#include <float.h>

namespace Spheral {

//------------------------------------------------------------------------------
// Access the artificial viscosity.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
ArtificialViscosityHandle<Dimension>&
GenericHydro<Dimension>::artificialViscosity() const {
  return mArtificialViscosity;
}

//------------------------------------------------------------------------------
// Access the CFL safety criteria.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
GenericHydro<Dimension>::cfl() const {
  return mCfl;
}

template<typename Dimension>
inline
void
GenericHydro<Dimension>::
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
GenericHydro<Dimension>::useVelocityMagnitudeForDt() const {
  return mUseVelocityMagnitudeForDt;
}

template<typename Dimension>
inline
void
GenericHydro<Dimension>::
useVelocityMagnitudeForDt(bool x) {
  mUseVelocityMagnitudeForDt = x;
}

//------------------------------------------------------------------------------
// Return the master neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GenericHydro<Dimension>::minMasterNeighbor() const {
  return mMinMasterNeighbor;
}

template<typename Dimension>
inline
int
GenericHydro<Dimension>::maxMasterNeighbor() const {
  return mMaxMasterNeighbor;
}

template<typename Dimension>
inline
double
GenericHydro<Dimension>::averageMasterNeighbor() const {
  return double(mSumMasterNeighbor)/(mNormMasterNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Return the coarse neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GenericHydro<Dimension>::minCoarseNeighbor() const {
  return mMinCoarseNeighbor;
}

template<typename Dimension>
inline
int
GenericHydro<Dimension>::maxCoarseNeighbor() const {
  return mMaxCoarseNeighbor;
}

template<typename Dimension>
inline
double
GenericHydro<Dimension>::averageCoarseNeighbor() const {
  return double(mSumCoarseNeighbor)/(mNormCoarseNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Return the refine neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GenericHydro<Dimension>::minRefineNeighbor() const {
  return mMinRefineNeighbor;
}

template<typename Dimension>
inline
int
GenericHydro<Dimension>::maxRefineNeighbor() const {
  return mMaxRefineNeighbor;
}

template<typename Dimension>
inline
double
GenericHydro<Dimension>::averageRefineNeighbor() const {
  return double(mSumRefineNeighbor)/(mNormRefineNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Return the actual neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GenericHydro<Dimension>::minActualNeighbor() const {
  return mMinActualNeighbor;
}

template<typename Dimension>
inline
int
GenericHydro<Dimension>::maxActualNeighbor() const {
  return mMaxActualNeighbor;
}

template<typename Dimension>
inline
double
GenericHydro<Dimension>::averageActualNeighbor() const {
  return double(mSumActualNeighbor)/(mNormActualNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Stored info about last timestep selection
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t
GenericHydro<Dimension>::
DTrank() const {
  return mDTrank;
}

template<typename Dimension>
inline
size_t
GenericHydro<Dimension>::
DTNodeList() const {
  return mDTNodeList;
}

template<typename Dimension>
inline
size_t
GenericHydro<Dimension>::
DTnode() const {
  return mDTnode;
}

template<typename Dimension>
inline
std::string
GenericHydro<Dimension>::
DTreason() const {
  return mDTreason;
}

//------------------------------------------------------------------------------
// Update the master neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
GenericHydro<Dimension>::
updateMasterNeighborStats(int numMaster) const {
  if (numMaster > 0) {
    mMinMasterNeighbor = std::min(mMinMasterNeighbor, numMaster);
    mMaxMasterNeighbor = std::max(mMaxMasterNeighbor, numMaster);
    mSumMasterNeighbor += numMaster;
    mNormMasterNeighbor++;
  }
}

//------------------------------------------------------------------------------
// Update the coarse neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
GenericHydro<Dimension>::
updateCoarseNeighborStats(int numNeighbor) const {
  if (numNeighbor > 0) {
    mMinCoarseNeighbor = std::min(mMinCoarseNeighbor, numNeighbor);
    mMaxCoarseNeighbor = std::max(mMaxCoarseNeighbor, numNeighbor);
    mSumCoarseNeighbor += numNeighbor;
    mNormCoarseNeighbor++;
  }
}

//------------------------------------------------------------------------------
// Update the refine neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
GenericHydro<Dimension>::
updateRefineNeighborStats(int numNeighbor) const {
  if (numNeighbor > 0) {
    mMinRefineNeighbor = std::min(mMinRefineNeighbor, numNeighbor);
    mMaxRefineNeighbor = std::max(mMaxRefineNeighbor, numNeighbor);
    mSumRefineNeighbor += numNeighbor;
    mNormRefineNeighbor++;
  }
}

//------------------------------------------------------------------------------
// Update the actual neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
GenericHydro<Dimension>::
updateActualNeighborStats(int numNeighbor) const {
  if (numNeighbor > 0) {
    mMinActualNeighbor = std::min(mMinActualNeighbor, numNeighbor);
    mMaxActualNeighbor = std::max(mMaxActualNeighbor, numNeighbor);
    mSumActualNeighbor += numNeighbor;
    mNormActualNeighbor++;
  }
}

}
