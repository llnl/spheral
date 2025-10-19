//---------------------------------Spheral++----------------------------------//
// GSPHFieldNames -- A collection of Field names specialized for GSPH module
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#ifndef _Spheral_GSPHFieldNames_
#define _Spheral_GSPHFieldNames_

#include <string>

namespace Spheral {

struct GSPHFieldNames {
  static const std::string nodalVelocity;
  static const std::string momentum;
  static const std::string thermalEnergy;
  static const std::string densityGradient;
  static const std::string specificThermalEnergyGradient;
  static const std::string pressureGradient;
  static const std::string deviatoricStressTensorGradient;
  static const std::string RiemannPressureGradient;
  static const std::string RiemannVelocityGradient;
  static const std::string RiemannDeviatoricStressTensorGradient;
  static const std::string pairMassFlux;
};

}

#endif
