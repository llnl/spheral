#-------------------------------------------------------------------------------
# PhysicsEvolvingMaterialLibrary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from SolidEquationOfState import *
from StrengthModel import *

@PYB11template("Dimension")
class PhysicsEvolvingMaterialLibrary(Physics,
                                     SolidEquationOfState,
                                     StrengthModel):
    """PhysicsEvolvingMaterialLibrary -- A base class to combine the Physics,
 EquationOfState, and StrengthModel for more sophisticated material
 libraries.

This class implements three distinct Spheral interfaces:
  EquationOfState : provides the ordinary Spheral EOS interface
  StrengthModel   : also answers the Strength questions about shear modulus 
                    and yield strength
  Physics         : Geodyn also needs to advance it's own internal state, so
                    this class implements the Physics interface to support
                    that."""

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ThirdRankTensor = typename %(Dimension)s::ThirdRankTensor;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
    using ScalarField = Field<%(Dimension)s, Scalar>;
"""


    def pyinit(self,
               referenceDensity = "const double",
               etamin = "const double",
               etamax = "const double",
               constants = "const PhysicalConstants&",
               minimumPressure = ("const double", "std::numeric_limits<double>::lowest()"),
               maximumPressure = ("const double", "std::numeric_limits<double>::max()"),
               minimumPressureDamage = ("const double", "0.0"),
               minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor"),
               externalPressure = ("const double", "0.0")):
        "Constructor"


