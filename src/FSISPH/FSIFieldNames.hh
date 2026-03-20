//---------------------------------Spheral++----------------------------------//
// FSIFieldNames -- Fields names used in FSI package
//----------------------------------------------------------------------------//
#ifndef _Spheral_FSIFieldNames_
#define _Spheral_FSIFieldNames_

#include <string>

namespace Spheral {

struct FSIFieldNames {
  const inline static std::string pressureGradient = "pressureGradient";
  const inline static std::string specificThermalEnergyGradient = "specificThermalEnergyGradient";
  const inline static std::string interfaceFlags = "interfaceFlags";
  const inline static std::string interfaceAreaVectors = "interfaceAreaVectors";
  const inline static std::string interfaceNormals = "interfaceNormals";
  const inline static std::string interfaceAngles = "interfaceAngles";
  const inline static std::string interfaceFraction = "interfaceFraction";
  const inline static std::string interfaceSmoothness = "interfaceSmoothness";
  const inline static std::string smoothedInterfaceNormals = "smoothedInterfaceNormals";
  const inline static std::string interfaceSmoothnessNormalization = "interfaceSmoothnessNormalization";
};

}

#endif
