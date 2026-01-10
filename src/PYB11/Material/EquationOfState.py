#-------------------------------------------------------------------------------
# EquationOfState abstract class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from EOSAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralMaterial")
class EquationOfState:

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef Field<%(Dimension)s, Scalar> ScalarField;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               constants = "const PhysicalConstants&",
               minimumPressure = "const double",
               maximumPressure = "const double",
               minPressureType = "const MaterialPressureMinType",
               externalPressure = ("const double", "0.0")):
        "EOS base constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def specificThermalEnergyForPressure(self,
                                         Ptarget = "const Scalar",
                                         rho = "const Scalar",
                                         epsMin = "const Scalar",
                                         epsMax = "const Scalar",
                                         epsTol = "const Scalar",
                                         Ptol = "const Scalar",
                                         maxIterations = ("const unsigned", "100"),
                                         verbose = ("const bool", "false")):
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def molecularWeight(self):
        "Optionally provide a molecular weight for an equation of state"
        return "Scalar"

    #...........................................................................
    # Methods
    @PYB11const
    def applyPressureLimits(self, P="const Scalar"):
        "Return the limited pressure"
        return "Scalar"

    #...........................................................................
    # Properties
    constants = PYB11property("const PhysicalConstants&", "constants", doc="The units choice")
    minimumPressure = PYB11property("double", "minimumPressure", "minimumPressure", doc="The minimum allowed pressure")
    maximumPressure = PYB11property("double", "maximumPressure", "maximumPressure", doc="The maximum allowed pressure")
    minimumPressureType = PYB11property("MaterialPressureMinType", "minimumPressureType", "minimumPressureType", doc="The algorithm for enforcing the minimum pressure")
    externalPressure = PYB11property("double", "externalPressure", "externalPressure")
    
#-------------------------------------------------------------------------------
# Add the virtual interface
#-------------------------------------------------------------------------------
PYB11inject(EOSAbstractMethods, EquationOfState, pure_virtual=True)
