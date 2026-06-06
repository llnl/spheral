#include "PhysicsEvolvingMaterialLibrary.hh"
#include "Field/Field.hh"
#include "Utilities/bisectSearch.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

#include <iostream>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PhysicsEvolvingMaterialLibrary<Dimension>::
PhysicsEvolvingMaterialLibrary(const double referenceDensity,
                               const double etamin,
                               const double etamax,
                               const PhysicalConstants& constants,
                               const double minimumPressure,
                               const double maximumPressure,
                               const double minimumPressureDamage,
                               const MaterialPressureMinType minPressureType,
                               const double externalPressure):
  Physics<Dimension>(),
  SolidEquationOfState<Dimension>(referenceDensity,
                                  etamin,
                                  etamax,
                                  constants,
                                  minimumPressure,
                                  maximumPressure,
                                  minimumPressureDamage,
                                  minPressureType,
                                  externalPressure),
  StrengthModel<Dimension>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PhysicsEvolvingMaterialLibrary<Dimension>::
~PhysicsEvolvingMaterialLibrary() {
}

}
