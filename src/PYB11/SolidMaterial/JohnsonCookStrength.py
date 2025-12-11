#-------------------------------------------------------------------------------
# JohnsonCookStrength
#-------------------------------------------------------------------------------
from PYB11Generator import *
from StrengthModel import *
from StrengthModelAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralSolidMaterial")
class JohnsonCookStrength(StrengthModel):
    "JohnsonCookStrength -- Implements the Johnson-Cook strength model."

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ScalarField = Field<%(Dimension)s, Scalar>;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               eos = "const SolidEquationOfState<%(Dimension)s>&",
               shearModulusModel = "const StrengthModel<%(Dimension)s>&",
               A = "const double",
               B = "const double",
               C = "const double",
               C4 = "const double",
               m = "const double",
               nhard = "const double",
               epsdot0 = "const double",
               epsdotmin = "const double",
               Tmelt = "const double",
               Troom = "const double",
               mu0 = "const double",
               shearModulusScaling = "const bool"):
        "Steinberg-Guinan strength model constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def providesSoundSpeed(self):
        return "bool"

    @PYB11virtual
    @PYB11const
    def soundSpeed(self,
                   soundSpeed = "Field<%(Dimension)s, Scalar>&",
                   density = "const Field<%(Dimension)s, Scalar>&",
                   specificThermalEnergy = "const Field<%(Dimension)s, Scalar>&",
                   pressure = "const Field<%(Dimension)s, Scalar>&",
                   fluidSoundSpeed = "const Field<%(Dimension)s, Scalar>&",
                   damage = "const Field<%(Dimension)s, SymTensor>&"):
        return "void"

    #...........................................................................
    # Properties
    A = PYB11property("double")
    B = PYB11property("double")
    C = PYB11property("double")
    C4 = PYB11property("double")
    m = PYB11property("double")
    nhard = PYB11property("double")
    epsdot0 = PYB11property("double")
    epsdotmin = PYB11property("double")
    Tmelt = PYB11property("double")
    Troom = PYB11property("double")
    mu0 = PYB11property("double")
    shearModulusScaling = PYB11property("bool")

#-------------------------------------------------------------------------------
# Inject Strength interface
#-------------------------------------------------------------------------------
PYB11inject(StrengthModelAbstractMethods, JohnsonCookStrength, virtual=True, pure_virtual=False)
