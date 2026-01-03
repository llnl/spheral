//---------------------------------Spheral++----------------------------------//
// SingularityEquationOfStateModels
//
// Constructors for different singularity EOS
//----------------------------------------------------------------------------//
#ifndef __Spheral_SingularityEquationOfStateModels_hh__
#define __Spheral_SingularityEquationOfStateModels_hh__

#include "SingularityEquationOfState.hh"

#include <singularity-eos/eos/eos.hpp>

namespace Spheral {

template<typename Dimension>
class SingularityGammaLawGas : public SingularityEquationOfState<Dimension> {
public:
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
                         const double externalPressure);
};

template<typename Dimension>
class SingularityExtendedVinet : public SingularityEquationOfState<Dimension> {
public:
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
                           const double externalPressure);
};

// Need singularity@1.9.0
// template<typename Dimension>
// class SingularityMieGruneisenLinear : public SingularityEquationOfState<Dimension> {
// public:
//   SingularityMieGruneisenLinear(const double rho0,
//                                 const double T0,
//                                 const double Cs,
//                                 const double s,
//                                 const double G0,
//                                 const double Cv0,
//                                 const double E0,
//                                 const double S0,
//                                 const PhysicalConstants& constants,
//                                 const SingularityLimitLevel level,
//                                 const double minimumPressure,
//                                 const double maximumPressure,
//                                 const double minimumPressureDamage,
//                                 const MaterialPressureMinType minPressureType,
//                                 const double externalPressure);
// };

template<typename Dimension>
class SingularitySpiner : public SingularityEquationOfState<Dimension> {
public:
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
                    const double externalPressure);
};

template<typename Dimension>
class SingularitySpinerT : public SingularityEquationOfState<Dimension> {
public:
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
                     const double externalPressure);
};

template<typename Dimension>
class SingularityEOSPAC : public SingularityEquationOfState<Dimension> {
public:
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
                    const double externalPressure);
};

}

#endif
