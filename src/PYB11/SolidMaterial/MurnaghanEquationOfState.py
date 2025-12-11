#-------------------------------------------------------------------------------
# MurnaghanEquationOfState
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SolidEquationOfState import *
from EOSAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralSolidMaterial")
class MurnaghanEquationOfState(SolidEquationOfState):
    """MurnaghanEquationOfState -- Murnaghan  equation of state.

  P(rho) = K/(n) * (eta^n - 1)
  eta = rho/rho0"""

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
               n = "const double",
               K = "const double",
               atomicWeight = "const double",
               constants = "const PhysicalConstants&",
               externalPressure = ("const double", "0.0"),
               minimumPressure = ("const double", "std::numeric_limits<double>::lowest()"),
               maximumPressure = ("const double", "std::numeric_limits<double>::max()"),
               minimumPressureDamage = ("const double", "0.0"),
               minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor")):
        "Murnaghan EOS"

    #...........................................................................
    # Methods
    @PYB11const
    def pressure(self,
                 massDensity = "const Scalar",
                 specificThermalEnergy = "const Scalar"):
        return "Scalar"

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
    n = PYB11property("double", "n", "n")
    K = PYB11property("double", "K", "K")

    atomicWeight = PYB11property("double", "atomicWeight", "atomicWeight")
    externalPressure = PYB11property("double", "externalPressure", "externalPressure")

#-------------------------------------------------------------------------------
# Inject EOS interface
#-------------------------------------------------------------------------------
PYB11inject(EOSAbstractMethods, MurnaghanEquationOfState, virtual=True)
