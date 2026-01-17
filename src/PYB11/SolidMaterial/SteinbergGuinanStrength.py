#-------------------------------------------------------------------------------
# SteinbergGuinanStrength
#-------------------------------------------------------------------------------
from PYB11Generator import *
from StrengthModel import *
from StrengthModelAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralSolidMaterial")
class SteinbergGuinanStrength(StrengthModel):
    """SteinbergGuinanStrength -- An implentation of StrengthModel returning constant
values for the shear modulus and yield strength."""

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ScalarField = Field<%(Dimension)s, Scalar>;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               eos = "const SolidEquationOfState<%(Dimension)s>&",
               G0 = "const double",
               Gmax = "const double",
               A = "const double",
               B = "const double",
               Y0 = "const double",
               Ymax = "const double",
               Yp = "const double",
               beta = "const double",
               gamma0 = "const double",
               nhard = "const double",
               coldEnergyFit = "const NinthOrderPolynomialFit&",
               meltEnergyFit = "const NinthOrderPolynomialFit&"):
        "Steinberg-Guinan strength model constructor"

    def pyinit1(self,
                eos = "const SolidEquationOfState<%(Dimension)s>&",
                G0 = "const double",
                A = "const double",
                B = "const double",
                Y0 = "const double",
                Ymax = "const double",
                Yp = "const double",
                beta = "const double",
                gamma0 = "const double",
                nhard = "const double",
                coldEnergyFit = "const NinthOrderPolynomialFit&",
                meltEnergyFit = "const NinthOrderPolynomialFit&"):
        "Steinberg-Guinan strength model constructor (without Gmax)"

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

    @PYB11virtual
    @PYB11const
    def meltSpecificEnergy(self,
                           meltSpecificEnergy = "Field<%(Dimension)s, Scalar>&",
                           density = "const Field<%(Dimension)s, Scalar>&",
                           specficThermalEnergy = "const Field<%(Dimension)s, Scalar>&"):
        return "void"

    @PYB11virtual
    @PYB11const
    def coldSpecificEnergy(self,
                           coldSpecificEnergy = "Field<%(Dimension)s, Scalar>&",
                           density = "const Field<%(Dimension)s, Scalar>&",
                           specficThermalEnergy = "const Field<%(Dimension)s, Scalar>&"):
        return "void"

    #...........................................................................
    # Methods
    @PYB11const
    def meltAttenuation(self,
                        rho = "const double",
                        eps = "const double"):
        "Melt attenuation."
        return "double"
    
    @PYB11const
    def computeTemperature(self,
                           temperature = "Field<%(Dimension)s, Scalar>&",
                           density = "const Field<%(Dimension)s, Scalar>&",
                           specificThermalEnergy = "const Field<%(Dimension)s, Scalar>&"):
        'Steinberg-Guinan "temperature".'
        return "void"

    #...........................................................................
    # Properties
    G0 = PYB11property("double")
    Gmax = PYB11property("double")
    A = PYB11property("double")
    B = PYB11property("double")
    Y0 = PYB11property("double")
    Ymax = PYB11property("double")
    Yp = PYB11property("double")
    beta = PYB11property("double")
    gamma0 = PYB11property("double")
    nhard = PYB11property("double")
    coldEnergyFit = PYB11property("const NinthOrderPolynomialFit&", returnpolicy="reference_internal")
    meltEnergyFit = PYB11property("const NinthOrderPolynomialFit&", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject Strength interface
#-------------------------------------------------------------------------------
PYB11inject(StrengthModelAbstractMethods, SteinbergGuinanStrength, virtual=True, pure_virtual=False)
