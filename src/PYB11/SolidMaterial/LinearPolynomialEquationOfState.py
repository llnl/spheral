#-------------------------------------------------------------------------------
# LinearPolynomialEquationOfState
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SolidEquationOfState import *
from EOSAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralSolidMaterial")
class LinearPolynomialEquationOfState(SolidEquationOfState):
    """LinearPolynomialEquationOfState -- An equation of state approximated by a
linear polynomial, i.e.:

  P(rho, e) = A0 + A1*mu + a2*mu^2 + a3*mu^3 + (B0 + B1*mu + B2*mu^2)*e
  mu = rho/rho0 - 1.0"""

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
               a0 = "const double",
               a1 = "const double",
               a2 = "const double",
               a3 = "const double",
               b0 = "const double",
               b1 = "const double",
               b2 = "const double",
               atomicWeight = "const double",
               constants = "const PhysicalConstants&",
               externalPressure = ("const double", "0.0"),
               minimumPressure = ("const double", "std::numeric_limits<double>::lowest()"),
               maximumPressure = ("const double", "std::numeric_limits<double>::max()"),
               minimumPressureDamage = ("const double", "0.0"),
               minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor")):
        "Linear-polynomial EOS"

    #...........................................................................
    # Methods
    @PYB11const
    @PYB11implementation("[](const LinearPolynomialEquationOfState<%(Dimension)s>& self, const Scalar& rho, const Scalar& eps) { return std::get<0>(self.pressureAndDerivs(rho, eps)); }")
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
    a0 = PYB11property("double", "a0", "a0")
    a1 = PYB11property("double", "a1", "a1")
    a2 = PYB11property("double", "a2", "a2")
    a3 = PYB11property("double", "a3", "a3")
    b0 = PYB11property("double", "b0", "b0")
    b1 = PYB11property("double", "b1", "b1")
    b2 = PYB11property("double", "b2", "b2")
    atomicWeight = PYB11property("double", "atomicWeight", "atomicWeight")
    externalPressure = PYB11property("double", "externalPressure", "externalPressure")

#-------------------------------------------------------------------------------
# Inject EOS interface
#-------------------------------------------------------------------------------
PYB11inject(EOSAbstractMethods, LinearPolynomialEquationOfState, virtual=True)
