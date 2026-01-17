#-------------------------------------------------------------------------------
# GruneisenEquationOfState
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SolidEquationOfState import *
from EOSAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralSolidMaterial")
class GruneisenEquationOfState(SolidEquationOfState):
    """GruneisenEquationOfState -- Gruneisen  equation of state.

Reference: Equation of State and Strength of Properties of Selected Materials
           Daniel J. Steinberg, UCRL-MA-106439, February 13, 1991"""

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef Field<%(Dimension)s, Scalar> ScalarField;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               referenceDensity = "const double",
               etamin = "const double",
               etamax = "const double",
               C0 = "const double",
               S1 = "const double",
               S2 = "const double",
               S3 = "const double",
               gamma0 = "const double",
               b = "const double",
               atomicWeight = "const double",
               constants = "const PhysicalConstants&",
               externalPressure = ("const double", "0.0"),
               minimumPressure = ("const double", "std::numeric_limits<double>::lowest()"),
               maximumPressure = ("const double", "std::numeric_limits<double>::max()"),
               minimumPressureDamage = ("const double", "0.0"),
               minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor")):
        "Gruneisen EOS"

    #...........................................................................
    # Methods
    @PYB11const
    @PYB11implementation("[](const GruneisenEquationOfState<%(Dimension)s>& self, const Scalar& rho, const Scalar& eps) { return std::get<0>(self.pressureAndDerivs(rho, eps)); }")
    def pressure(self,
                 massDensity = "const Scalar",
                 specificThermalEnergy = "const Scalar"):
        return "Scalar"

    @PYB11const
    def pressureAndDerivs(self,
                          massDensity = "const Scalar",
                          specificThermalEnergy = "const Scalar"):
        return "std::tuple<Scalar, Scalar, Scalar>"

    @PYB11const
    def temperature(self,
                    massDensity = "const Scalar",
                    specificThermalEnergy = "const Scalar"):
        return "Scalar"

    @PYB11const
    def specificThermalEnergy(self,
                              massDensity = "const Scalar",
                              temperature = "const Scalar"):
        return "Scalar"

    @PYB11const
    def specificHeat(self,
                     massDensity = "const Scalar",
                     temperature = "const Scalar"):
        return "Scalar"

    @PYB11const
    def soundSpeed(self,
                   massDensity = "const Scalar",
                   specificThermalEnergy = "const Scalar"):
        return "Scalar"

    @PYB11const
    def gamma(self,
              massDensity = "const Scalar",
              specificThermalEnergy = "const Scalar"):
        return "Scalar"

    @PYB11const
    def bulkModulus(self,
                    massDensity = "const Scalar",
                    specificThermalEnergy = "const Scalar"):
        return "Scalar"

    @PYB11const
    def entropy(self,
                massDensity = "const Scalar",
                specificThermalEnergy = "const Scalar"):
        return "Scalar"

    @PYB11const
    def computeDPDrho(self,
                      massDensity = "const Scalar",
                      specificThermalEnergy = "const Scalar"):
        "Compute the derivative of the pressure with respect to the density."
        return "double"

    #...........................................................................
    # Properties
    C0 = PYB11property("double", "C0", "C0")
    S1 = PYB11property("double", "S1", "S1")
    S2 = PYB11property("double", "S2", "S2")
    S3 = PYB11property("double", "S3", "S3")
    gamma0 = PYB11property("double", "gamma0", "gamma0")
    b = PYB11property("double", "b", "b")
    Cv = PYB11property("double", "Cv")
    atomicWeight = PYB11property("double", "atomicWeight", "atomicWeight")
    energyMultiplier = PYB11property("double", "energyMultiplier", "energyMultiplier",
                                     doc="""Option to scale the thermal energy term by.  This is mostly useful for test problems
where you want to make the Gruneisen independent of energy.""")
    externalPressure = PYB11property("double", "externalPressure", "externalPressure")

#-------------------------------------------------------------------------------
# Inject EOS interface
#-------------------------------------------------------------------------------
PYB11inject(EOSAbstractMethods, GruneisenEquationOfState, virtual=True)
