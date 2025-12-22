#-------------------------------------------------------------------------------
# TillotsonEquationOfState
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SolidEquationOfState import *
from EOSAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralSolidMaterial")
class TillotsonEquationOfState(SolidEquationOfState):
    """TillotsonEquationOfState -- Tillotson  equation of state.

This equation of state is designed to represent metallic materials
over a range of pressure and density -- spanning solid, liquid, and
vapor states.

Reference: Tillotson 1962
    """

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
               etamin_solid = "const double",
               etamax_solid = "const double",
               a = "const double",
               b = "const double",
               A = "const double",
               B = "const double",
               alpha = "const double",
               beta = "const double",
               eps0 = "const double",
               epsLiquid = "const double",
               epsVapor = "const double",
               atomicWeight = "const double",
               constants = "const PhysicalConstants&",
               externalPressure = ("const double", "0.0"),
               minimumPressure = ("const double", "std::numeric_limits<double>::lowest()"),
               maximumPressure = ("const double", "std::numeric_limits<double>::max()"),
               minimumPressureDamage = ("const double", "0.0"),
               minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor")):
        "Tillotson EOS"

    #...........................................................................
    # Methods
    # Methods
    @PYB11const
    @PYB11implementation("[](const TillotsonEquationOfState<%(Dimension)s>& self, const Scalar& rho, const Scalar& eps) { return std::get<0>(self.pressureAndDerivs(rho, eps)); }")
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

    @PYB11const
    def computePhi(self,
                   eta = "const double&",
                   eps = "const double&"):
        return "double"

    @PYB11const
    def computeP1(self,
                  mu = "const double&",
                  P2 = "const double&"):
        return "double"

    @PYB11const
    def computeP2(self,
                  phi = "const double&",
                  mu = "const double&",
                  rho = "const double&",
                  eps = "const double&"):
        return "double"

    @PYB11const
    def computeP4(self,
                  phi = "const double&",
                  mu = "const double&",
                  eta = "const double&",
                  rho = "const double&",
                  eps = "const double&"):
        return "double"

    @PYB11const
    def compute_dphidrho_eps(self,
                             rho0 = "const double&",
                             eta = "const double&",
                             eps = "const double&"):
        return "double"

    @PYB11const
    def compute_dP1drho_eps(self,
                            rho0 = "const double&",
                            mu = "const double&",
                            eps = "const double& dP2drho"):
        return "double"

    @PYB11const
    def compute_dP2drho_eps(self,
                            phi = "const double&",
                            dphidrho_eps = "const double&",
                            rho0 = "const double&",
                            rho = "const double&",
                            eps = "const double&"):
        return "double"

    @PYB11const
    def compute_dP4drho_eps(self,
                            phi = "const double&",
                            dphidrho_eps = "const double&",
                            rho0 = "const double&",
                            eta = "const double&",
                            mu = "const double&",
                            rho = "const double&",
                            eps = "const double&"):
        return "double"

    @PYB11const
    def compute_dphideps_rho(self,
                             eta = "const double&",
                             eps = "const double&"):
        return "double"

    @PYB11const
    def compute_dP2deps_rho(self,
                            phi = "const double&",
                            dphideps_rho = "const double&",
                            rho = "const double&",
                            eps = "const double&"):
        return "double"

    @PYB11const
    def compute_dP4deps_rho(self,
                            phi = "const double&",
                            dphideps_rho = "const double&",
                            eta = "const double&",
                            rho = "const double&",
                            eps = "const double&"):
        return "double"

    #...........................................................................
    # Properties
    etamin_solid = PYB11property("double", "etamin_solid", "etamin_solid")
    etamax_solid = PYB11property("double", "etamax_solid", "etamax_solid")
    a = PYB11property("double", "a", "a")
    b = PYB11property("double", "b", "b")
    A = PYB11property("double", "A", "A")
    B = PYB11property("double", "B", "B")
    alpha = PYB11property("double", "alpha", "alpha")
    beta = PYB11property("double", "beta", "beta")
    eps0 = PYB11property("double", "eps0", "eps0")
    epsLiquid = PYB11property("double", "epsLiquid", "epsLiquid")
    epsVapor = PYB11property("double", "epsVapor", "epsVapor")

    atomicWeight = PYB11property("double", "atomicWeight", "atomicWeight")

#-------------------------------------------------------------------------------
# Inject EOS interface
#-------------------------------------------------------------------------------
PYB11inject(EOSAbstractMethods, TillotsonEquationOfState, virtual=True)
