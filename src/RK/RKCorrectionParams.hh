//Parameters used for RK kernels. 
#ifndef __Spheral_RKCorrectionParams_hh__
#define __Spheral_RKCorrectionParams_hh__

#include "VoronoiCells/VolumeType.hh"

#include <iostream>

//Enumerated type for the corrected Kernels
namespace Spheral {

enum class RKOrder : int {//Used to assign the order of the corrections
  ZerothOrder = 0,
  LinearOrder = 1,
  QuadraticOrder = 2,
  CubicOrder = 3,
  QuarticOrder = 4,
  QuinticOrder = 5,
  SexticOrder = 6,
  SepticOrder = 7,
};

// Backward compatibility
using RKVolumeType = VolumeType;

//------------------------------------------------------------------------------
// ostream operator for RKOrder
//------------------------------------------------------------------------------
inline
std::ostream& operator<<(std::ostream& os, const RKOrder& x) {
  switch(x) {
  case RKOrder::ZerothOrder:    os << "ZerothOrder"; break;
  case RKOrder::LinearOrder:    os << "LinearOrder"; break;
  case RKOrder::QuadraticOrder: os << "QuadraticOrder"; break;
  case RKOrder::CubicOrder:     os << "CubicOrder"; break;
  case RKOrder::QuarticOrder:   os << "QuarticOrder"; break;
  case RKOrder::QuinticOrder:   os << "QuinticOrder"; break;
  case RKOrder::SexticOrder:    os << "SexticOrder"; break;
  case RKOrder::SepticOrder:    os << "SepticOrder"; break;
  default:                      os << "Unknown RKOrder";
  }
  return os;
}

}
#endif
