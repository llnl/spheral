//---------------------------------Spheral++----------------------------------//
// PhysicsEvolvingMaterialLibrary -- A base class to combine the Physics,
//  EquationOfState, and StrengthModel for more sophisticated material
//  libraries.
//
// This class implements three distinct Spheral interfaces:
//   EquationOfState : provides the ordinary Spheral EOS interface
//   StrengthModel   : also answers the Strength questions about shear modulus 
//                     and yield strength
//   Physics         : Geodyn also needs to advance it's own internal state, so
//                     this class implements the Physics interface to support
//                     that.
//
// Created by JMO, Thu Mar 31 16:11:26 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_PhysicsEvolvingMaterialLibrary_hh__
#define __Spheral_PhysicsEvolvingMaterialLibrary_hh__

#include "SolidMaterial/SolidEquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "Physics/Physics.hh"

namespace Spheral {

template<typename Dimension>
class PhysicsEvolvingMaterialLibrary: 
    public Physics<Dimension>,
    public SolidEquationOfState<Dimension>,
    public StrengthModel<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors.
  PhysicsEvolvingMaterialLibrary(const double referenceDensity,
                                 const double etamin,
                                 const double etamax,
                                 const PhysicalConstants& constants,
                                 const double minimumPressure,
                                 const double maximumPressure,
                                 const double minimumPressureDamage,
                                 const MaterialPressureMinType minPressureType,
                                 const double externalPressure);
  virtual ~PhysicsEvolvingMaterialLibrary();

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#endif
