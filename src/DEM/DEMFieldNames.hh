//---------------------------------Spheral++----------------------------------//
// DEMFieldNames -- A collection of standard Field names for the DEM 
// physics package.
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//
#ifndef _Spheral_DEMFieldNames_
#define _Spheral_DEMFieldNames_

#include <string>

namespace Spheral {

struct DEMFieldNames {
  static inline const std::string momentOfInertia = "moment of inertia";
  static inline const std::string particleRadius = "particle radius";
  static inline const std::string compositeParticleIndex = "composite particle flags";
  static inline const std::string angularVelocity = "angular velocity";
  static inline const std::string uniqueIndices = "unique indices";
  static inline const std::string isActiveContact = "bool indentifying active contacts";
  static inline const std::string neighborIndices = "unique neighbor indices";
  static inline const std::string shearDisplacement = "shear displacement";
  static inline const std::string rollingDisplacement = "rolling displacement";
  static inline const std::string torsionalDisplacement = "torsional displacement";
  static inline const std::string equilibriumOverlap = "equilibrium overlap";
  static inline const std::string maximumOverlap = "maximum overlap";
  static inline const std::string solidBoundaries = "solid boundaries";
  static inline const std::string solidBoundaryPolicy = "solid boundary policy";
};

}

#endif
