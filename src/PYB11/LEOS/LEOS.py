#-------------------------------------------------------------------------------
# LEOS
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SolidEquationOfState import *
from EOSAbstractMethods import *

@PYB11template("Dimension")
class LEOS(SolidEquationOfState):
    """LEOS

The Livermore equation of state.
"""

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using ScalarField = Field<%(Dimension)s, Scalar>;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               materialNumber = "const int",
               constants = "const PhysicalConstants&",
               externalPressure = ("const double", "0.0"),
               minimumPressure = ("const double", "std::numeric_limits<double>::lowest()"),
               maximumPressure = ("const double",  "std::numeric_limits<double>::max()"),
               minimumPressureDamage = ("const double", "0.0"),
               minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor"),
               dbname = ("const std::string", '"leos"'),
               leosFileFormat = ("const std::string", '""'),
               atomicWeight = ("const double", "0.0")):
        "LEOS constructor"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def setMeltTemperature(self,
                           meltTemperature = "ScalarField&",
                           massDensity = "const ScalarField&",
                           specificThermalEnergy = "const ScalarField&"):
        return "void"

    @PYB11const
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
    @PYB11pycppname("gamma")
    def gamma1(self,
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
    def meltTemperature(self,
                        massDensity = "const Scalar",
                        specificThermalEnergy = "const Scalar"):
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def specificThermalEnergyForPressure(self,
                                         Ptarget = "const Scalar",
                                         rho = "const Scalar",
                                         epsMin = ("const Scalar", "0.0"),
                                         epsMax = ("const Scalar", "0.0"),
                                         epsTol = ("const Scalar", "0.0"),
                                         Ptol = ("const Scalar", "0.0"),
                                         maxIterations = ("const unsigned", "100"),
                                         verbose = ("const bool", "false")):
        "Override the default EOS method of searching for the specific energy corresponding to a pressure with an LEOS inverse lookup"
        return "Scalar"

    #...........................................................................
    # Properties
    materialNumber = PYB11property()
    databaseName = PYB11property()
    referenceTemperature = PYB11property()
    atomicWeight = PYB11property()
    K0 = PYB11property()
    descriptor = PYB11property()
    
#-------------------------------------------------------------------------------
# Add the virtual interface
#-------------------------------------------------------------------------------
PYB11inject(EOSAbstractMethods, LEOS, virtual=True)
