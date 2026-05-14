//------------------------------------------------------------------------------
// SPHRZUtilities
//
// Common methods to SPHRZ and SolidSPHRZ.
//------------------------------------------------------------------------------
#ifndef __Spheral_SPHRZUtilities__
#define __Spheral_SPHRZUtilities__

#include "Physics/GenericHydro.hh"
#include "SPH/computeSPHSumMassDensity.hh"
#include "SPH/correctSPHSumMassDensity.hh"
#include "SPH/computeSumVoronoiCellMassDensity.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/AxisymmetricMassDensityPolicy.hh"
#include "Hydro/RZNonSymmetricSpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"

#include <limits>

namespace Spheral {
namespace SPHRZUtilities {

//------------------------------------------------------------------------------
// Update the real mass density based on the areal (RZ) density (Fields)
//------------------------------------------------------------------------------
inline
void
computeMassDensityFromArealMassDensity(Field<Dim<2>, Dim<2>::Scalar>& massDensity,
                                       const Field<Dim<2>, Dim<2>::Scalar>& massDensityRZ,
                                       const Field<Dim<2>, Dim<2>::Scalar>& mass,
                                       const Field<Dim<2>, Dim<2>::Scalar>& massRZ,
                                       const Field<Dim<2>, Dim<2>::Vector>& position,
                                       const Field<Dim<2>, Dim<2>::SymTensor>& H,
                                       const Dim<2>::Scalar rhoMin = std::numeric_limits<Dim<2>::Scalar>::min(),
                                       const Dim<2>::Scalar rhoMax = std::numeric_limits<Dim<2>::Scalar>::max()) {
  const auto n = massDensity.numInternalElements();
  for (auto i = 0u; i < n; ++i) {
    CHECK(massDensityRZ(i) > 0.0);
    const auto& posi = position(i);
    const auto& Hi = H(i);
    const auto  ri = std::abs(posi.y());
    const auto  zetai = std::abs((Hi*posi).y());
    const auto  hri = ri*safeInv(zetai);
    CHECK(hri >= 0.0);
    const auto Ai = massRZ(i)/massDensityRZ(i);
    const auto Vi = 2.0*M_PI*ri*Ai;
    const auto Vmin = 0.001*FastMath::pow3(hri);
    const auto rho3D = std::clamp(mass(i)*safeInvVar(Vi, Vmin), rhoMin, rhoMax);
    const auto w0 = std::clamp(ri/(0.1*hri) - 0.05, 0.0, 1.0);
    massDensity(i) = w0*rho3D + (1.0 - w0)*massDensityRZ(i);
    // CHECK(Vi > 0.0);
    // CHECK(Vmin > 0.0);
    // massDensity(i) = std::clamp(mass(i)*safeInvVar(Vi, Vmin), rhoMin, rhoMax);
  }
}
  
//------------------------------------------------------------------------------
// Update the real mass density based on the areal (RZ) density (FieldLists)
//------------------------------------------------------------------------------
inline
void
computeMassDensityFromArealMassDensity(FieldList<Dim<2>, Dim<2>::Scalar>& massDensity,
                                       const FieldList<Dim<2>, Dim<2>::Scalar>& massDensityRZ,
                                       const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                       const FieldList<Dim<2>, Dim<2>::Scalar>& massRZ,
                                       const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                       const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                                       const Dim<2>::Scalar rhoMin = std::numeric_limits<Dim<2>::Scalar>::min(),
                                       const Dim<2>::Scalar rhoMax = std::numeric_limits<Dim<2>::Scalar>::max()) {
  for (auto k = 0u; k < massDensity.numFields(); ++k) {
    computeMassDensityFromArealMassDensity(*massDensity[k],
                                           *massDensityRZ[k],
                                           *mass[k],
                                           *massRZ[k],
                                           *position[k],
                                           *H[k],
                                           rhoMin,
                                           rhoMax);
  }
}
  
//------------------------------------------------------------------------------
// preStepInitialize
//------------------------------------------------------------------------------
template<typename HydroPackage>
inline
void
preStepInitialize(const HydroPackage& hydro,
                  const DataBase<Dim<2>>& dataBase, 
                  State<Dim<2>>& state,
                  StateDerivatives<Dim<2>>& derivs) {

  using Vector = Dim<2>::Vector;
  using SymTensor = Dim<2>::SymTensor;

  // Grab the state we need
  const auto& WT = hydro.kernel();
  const auto  densityUpdate = hydro.densityUpdate();
  const auto  sumMassDensityOverAllNodeLists = hydro.sumMassDensityOverAllNodeLists();
  const auto& boundaries = hydro.boundaryConditions();
  const auto& connectivityMap = state.connectivityMap();
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto  massRZ = state.fields(HydroFieldNames::massRZ, 0.0);
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero());
  auto        massDensityRZ = state.fields(HydroFieldNames::massDensityRZ, 0.0);
  auto        massDensity = state.fields(HydroFieldNames::massDensity, 0.0);

  switch(densityUpdate) {

  case MassDensityType::IntegrateDensity:
    break;

  case MassDensityType::RigorousSumDensity:
  case MassDensityType::CorrectedSumDensity:
    {
      computeSPHSumMassDensity(connectivityMap, WT, sumMassDensityOverAllNodeLists, position, massRZ, H, massDensityRZ);
      if (densityUpdate == MassDensityType::CorrectedSumDensity) {
        for (auto* boundPtr: boundaries) boundPtr->applyFieldListGhostBoundary(massDensityRZ);
        for (auto* boundPtr: boundaries) boundPtr->finalizeGhostBoundary();
        correctSPHSumMassDensity(connectivityMap, WT, sumMassDensityOverAllNodeLists, position, massRZ, H, massDensityRZ);
      }
    }
    break;

  case MassDensityType::VoronoiCellDensity:
    {
      hydro.updateVolume(state, false);
      const auto volume = state.fields(HydroFieldNames::volume, 0.0);
      massDensityRZ = massRZ / volume;
    }
    break;

  case MassDensityType::SumVoronoiCellDensity:
    {
      hydro.updateVolume(state, true);
      const auto volume = state.fields(HydroFieldNames::volume, 0.0);
      computeSumVoronoiCellMassDensity(connectivityMap, WT, position, massRZ, volume, H, massDensityRZ);
    }
    break;

  default:
    VERIFY2(false, "preStepInitialize did not handle a density update choice : " << static_cast<int>(densityUpdate));
    break;
  }

  // Update the real mass density based on the areal (RZ) density
  computeMassDensityFromArealMassDensity(massDensity, massDensityRZ, mass, massRZ, position, H);
  for (auto* boundPtr: boundaries) {
    boundPtr->applyFieldListGhostBoundary(massDensityRZ);
    boundPtr->applyFieldListGhostBoundary(massDensity);
  }
  for (auto* boundPtr: boundaries) boundPtr->finalizeGhostBoundary();
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename HydroPackage>
inline
void
registerState(const HydroPackage& hydro,
              DataBase<Dim<2>>& dataBase,
              State<Dim<2>>& state,
              FieldList<Dim<2>, Dim<2>::Scalar>& massRZ,
              FieldList<Dim<2>, Dim<2>::Scalar>& massDensityRZ) {

  using Dimension = Dim<2>;
  using Scalar = Dimension::Scalar;

  // RZ mass
  state.enroll(massRZ);

  // Set the mass density policies
  auto rho = state.fields(HydroFieldNames::massDensity, 0.0);
  for (auto [nodeListi, fluidNodeListPtr]: enumerate(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
    state.enroll(*massDensityRZ[nodeListi], make_policy<IncrementBoundedState<Dimension, Scalar>>(fluidNodeListPtr->rhoMin(),
                                                                                                  fluidNodeListPtr->rhoMax()));
    state.enroll(*rho[nodeListi], make_policy<AxisymmetricMassDensityPolicy>(fluidNodeListPtr->rhoMin(),
                                                                             fluidNodeListPtr->rhoMax()));
  }

  // // Reregister the plastic strain policy to the RZ specialized version
  // // that accounts for the theta-theta component of the stress.  Also the deviatoric stress.
  // auto ps = state.fields(SolidFieldNames::plasticStrain, 0.0);
  // auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero());
  // PolicyPointer plasticStrainPolicy(new RZPlasticStrainPolicy());
  // PolicyPointer deviatoricStressPolicy(new DeviatoricStressPolicy<Dimension>(false));
  // state.enroll(ps, plasticStrainPolicy);
  // state.enroll(S, deviatoricStressPolicy);

  // We have to choose either compatible or total energy evolution.
  const auto compatibleEnergy = hydro.compatibleEnergyEvolution();
  const auto evolveTotalEnergy = hydro.evolveTotalEnergy();
  VERIFY2(not (compatibleEnergy and evolveTotalEnergy),
          "SPH error : you cannot simultaneously use both compatibleEnergyEvolution and evolveTotalEnergy");

  // Register the specific thermal energy.
  // Note in RZ we require the specific thermal energy go before the position so we can use the r position
  // during update.  This is why we make position update dependent on the thermal energy in SPHBase.
  auto specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  if (compatibleEnergy) {
    state.enroll(specificThermalEnergy, make_policy<RZNonSymmetricSpecificThermalEnergyPolicy>(dataBase));

  } else if (evolveTotalEnergy) {
    // If we're doing total energy, we register the specific energy to advance with the
    // total energy policy.
    state.enroll(specificThermalEnergy, make_policy<SpecificFromTotalThermalEnergyPolicy<Dimension>>());

  } else {
    // Otherwise we're just time-evolving the specific energy.
    state.enroll(specificThermalEnergy, make_policy<IncrementState<Dimension, Scalar>>());
  }
}

}
}

#endif
