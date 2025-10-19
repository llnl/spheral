//---------------------------------Spheral++----------------------------------//
// GSPHFieldNames -- A collection of Field names specialized for GSPH module
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#include "GSPHFieldNames.hh"

const std::string Spheral::GSPHFieldNames::nodalVelocity = "velocity of node";
const std::string Spheral::GSPHFieldNames::momentum = "momentum";
const std::string Spheral::GSPHFieldNames::thermalEnergy = "thermal energy";
const std::string Spheral::GSPHFieldNames::densityGradient = "density gradient";
const std::string Spheral::GSPHFieldNames::specificThermalEnergyGradient = "specific thermal energy gradient";
const std::string Spheral::GSPHFieldNames::pressureGradient = "pressure gradient";
const std::string Spheral::GSPHFieldNames::deviatoricStressTensorGradient = "deviatoric stress tensor gradient";
const std::string Spheral::GSPHFieldNames::RiemannPressureGradient = "Riemann solvers pressure gradient";
const std::string Spheral::GSPHFieldNames::RiemannVelocityGradient = "Riemann solvers velocity gradient";
const std::string Spheral::GSPHFieldNames::RiemannDeviatoricStressTensorGradient = "Riemann solvers deviatoric stress tensor gradient";
const std::string Spheral::GSPHFieldNames::pairMassFlux = "pairwise mass flux";