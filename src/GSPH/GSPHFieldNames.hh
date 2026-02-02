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
  const inline static std::string nodalVelocity = "velocity of node";
  const inline static std::string momentum = "momentum";
  const inline static std::string thermalEnergy = "thermal energy";
  const inline static std::string densityGradient = "density gradient";
  const inline static std::string pressureGradient = "pressure gradient";
  const inline static std::string deviatoricStressTensorGradient = "deviatoric stress tensor gradient";
  const inline static std::string RiemannPressureGradient = "Riemann solvers pressure gradient";
  const inline static std::string RiemannVelocityGradient = "Riemann solvers velocity gradient";
  const inline static std::string RiemannDeviatoricStressTensorGradient = "Riemann solvers deviatoric stress tensor gradient";
  const inline static std::string pairMassFlux = "pairwise mass flux";
};

}

#endif
