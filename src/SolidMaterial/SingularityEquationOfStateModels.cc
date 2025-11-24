//---------------------------------Spheral++----------------------------------//
// SingularityEquationOfStateModels
//
// Constructors for different singularity EOS
//----------------------------------------------------------------------------//
#include "SingularityEquationOfStateModels.hh"
#include "Material/PhysicalConstants.hh"
#include "Utilities/SpheralFunctions.hh"

namespace Spheral {

namespace { // anonymous

template<typename EOS, typename... Args>
SingularityEOSType newEOS(const PhysicalConstants& constants,
                          Args&&... args) {
  return EOS(std::forward<Args>(args)...);
  // return singularity::UnitSystem<EOS>(
  //   EOS(std::forward<Args>(args)...),
  //   singularity::eos_units_init::LengthTimeUnitsInit(),
  //   constants.unitTimeSec(),
  //   1000 * constants.unitMassKg(),
  //   100 * constants.unitLengthMeters(),
  //   constants.unitTemperatureKelvin());
}

double getCvCGS(const double gamma,
                const double mu) {
  // Get cv in CGS units
  PhysicalConstants constants(0.01, 0.001, 1.0, 1.0);
  const auto kB = constants.kB();
  const auto mp = constants.protonMass();
  return kB / ((gamma - 1) * mu * mp);
}

} // namespace anonymous

template<typename Dimension>
SingularityGammaLawGas<Dimension>::
SingularityGammaLawGas(const double gamma,
                       const double mu,
                       const PhysicalConstants& constants,
                       const SingularityLimitLevel level,
                       const double referenceDensity,
                       const double etamin,
                       const double etamax,
                       const double minimumPressure,
                       const double maximumPressure,
                       const double minimumPressureDamage,
                       const MaterialPressureMinType minPressureType,
                       const double externalPressure):
  SingularityEquationOfState<Dimension>(
    newEOS<singularity::IdealGas>(constants, gamma - 1, getCvCGS(gamma, mu)),
    constants, level, referenceDensity, etamin, etamax,
    minimumPressure, maximumPressure, minimumPressureDamage,
    minPressureType, externalPressure) {
}

template<typename Dimension>
SingularityExtendedVinet<Dimension>::
SingularityExtendedVinet(const double rho0,
                         const double T0,
                         const double B0,
                         const double BP0,
                         const double A0,
                         const double Cv0,
                         const double E0,
                         const double S0,
                         const std::vector<double>& coeffs,
                         const PhysicalConstants& constants,
                         const SingularityLimitLevel level,
                         const double referenceDensity,
                         const double etamin,
                         const double etamax,
                         const double minimumPressure,
                         const double maximumPressure,
                         const double minimumPressureDamage,
                         const MaterialPressureMinType minPressureType,
                         const double externalPressure):
  SingularityEquationOfState<Dimension>(
    newEOS<singularity::Vinet>(constants, rho0, T0, B0, BP0, A0, Cv0, E0, S0, &coeffs[0]),
    constants, level, referenceDensity, etamin, etamax,
    minimumPressure, maximumPressure, minimumPressureDamage,
    minPressureType, externalPressure) {
}

// template<typename Dimension>
// SingularityMieGruneisenLinear<Dimension>::
// SingularityMieGruneisenLinear(const double rho0,
//                               const double T0,
//                               const double Cs,
//                               const double s,
//                               const double G0,
//                               const double Cv0,
//                               const double E0,
//                               const double S0,
//                               const PhysicalConstants& constants,
//                               const SingularityLimitLevel level,
//                               const double minimumPressure,
//                               const double referenceDensity,
//                               const double etamin,
//                               const double etamax,
//                               const double maximumPressure,
//                               const double minimumPressureDamage,
//                               const MaterialPressureMinType minPressureType,
//                               const double externalPressure):
//   SingularityEquationOfState<Dimension>(
//     newEOS<singularity::MGUsup>(constants, rho0, T0, Cs, s, G0, Cv0, E0, S0),
//     constants, level, referenceDensity, etamin, etamax,
//     minimumPressure, maximumPressure, minimumPressureDamage,
//     minPressureType, externalPressure) {
// }

template<typename Dimension>
SingularitySpiner<Dimension>::
SingularitySpiner(const std::string &filename,
                  int matid,
                  const PhysicalConstants& constants,
                  const SingularityLimitLevel level,
                  const double referenceDensity,
                  const double etamin,
                  const double etamax,
                  const double minimumPressure,
                  const double maximumPressure,
                  const double minimumPressureDamage,
                  const MaterialPressureMinType minPressureType,
                  const double externalPressure):
  SingularityEquationOfState<Dimension>(
    newEOS<singularity::SpinerEOSDependsRhoSie>(constants, filename, matid),
    constants, level, referenceDensity, etamin, etamax,
    minimumPressure, maximumPressure, minimumPressureDamage,
    minPressureType, externalPressure) {
}

template<typename Dimension>
SingularitySpinerT<Dimension>::
SingularitySpinerT(const std::string &filename,
                   int matid,
                   const PhysicalConstants& constants,
                   const SingularityLimitLevel level,
                   const double referenceDensity,
                   const double etamin,
                   const double etamax,
                   const double minimumPressure,
                   const double maximumPressure,
                   const double minimumPressureDamage,
                   const MaterialPressureMinType minPressureType,
                   const double externalPressure):
  SingularityEquationOfState<Dimension>(
    newEOS<singularity::SpinerEOSDependsRhoT>(constants, filename, matid),
    constants, level, referenceDensity, etamin, etamax,
    minimumPressure, maximumPressure, minimumPressureDamage,
    minPressureType, externalPressure) {
}

template<typename Dimension>
SingularityEOSPAC<Dimension>::
SingularityEOSPAC(int matid,
                  const PhysicalConstants& constants,
                  const SingularityLimitLevel level,
                  const double referenceDensity,
                  const double etamin,
                  const double etamax,
                  const double minimumPressure,
                  const double maximumPressure,
                  const double minimumPressureDamage,
                  const MaterialPressureMinType minPressureType,
                  const double externalPressure):
  SingularityEquationOfState<Dimension>(
    newEOS<singularity::EOSPAC>(constants, matid, true), // invert at setup
    constants, level, referenceDensity, etamin, etamax,
    minimumPressure, maximumPressure, minimumPressureDamage,
    minPressureType, externalPressure) {
}

} // namespace Spheral
