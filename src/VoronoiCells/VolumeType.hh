//---------------------------------Spheral++----------------------------------//
// VolumeType
//
// Enumeration for volume computation methods.
//----------------------------------------------------------------------------//
#ifndef __Spheral_VolumeType_hh__
#define __Spheral_VolumeType_hh__

namespace Spheral {

enum class VolumeType : int {
  MassOverDensity = 0,
  SumVolume       = 1,
  VoronoiVolume   = 2,
  HullVolume      = 3,
  HVolume         = 4,
  // Backward-compatible aliases
  RKMassOverDensity = MassOverDensity,
  RKSumVolume       = SumVolume,
  RKVoronoiVolume   = VoronoiVolume,
  RKHullVolume      = HullVolume,
};

} // end namespace Spheral

#endif
