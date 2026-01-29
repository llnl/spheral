#-------------------------------------------------------------------------------
# OsborneEquationOfState
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SolidEquationOfState import *
from EOSAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralSolidMaterial")
class OsborneEquationOfState(SolidEquationOfState):
    """OsborneEquationOfState -- Osborne  equation of state.
Reference: PAGOSA Physics manual, LA-14425-M"""

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
               a1 = "const double",
               a2pos = "const double",
               a2neg = "const double",
               b0 = "const double",
               b1 = "const double",
               b2pos = "const double",
               b2neg = "const double",
               c0 = "const double",
               c1 = "const double",
               c2pos = "const double",
               c2neg = "const double",
               E0 = "const double",
               atomicWeight = "const double",
               constants = "const PhysicalConstants&",
               externalPressure = ("const double", "0.0"),
               minimumPressure = ("const double", "std::numeric_limits<double>::lowest()"),
               maximumPressure = ("const double", "std::numeric_limits<double>::max()"),
               minimumPressureDamage = ("const double", "0.0"),
               minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor")):
        "Osborne EOS"

    #...........................................................................
    # Methods
    @PYB11const
    def DPDrho(self,
               massDensity = "const double",
               specificThermalEnergy = "const double"):
        "Compute the derivative of the pressure with respect to the density."
        return "double"

    #...........................................................................
    # Properties
    a1 = PYB11property("double", "a1", "a1")
    a2pos = PYB11property("double", "a2pos", "a2pos")
    a2neg = PYB11property("double", "a2neg", "a2neg")
    b0 = PYB11property("double", "b0", "b0")
    b1 = PYB11property("double", "b1", "b1")
    b2pos = PYB11property("double", "b2pos", "b2pos")
    b2neg = PYB11property("double", "b2neg", "b2neg")
    c0 = PYB11property("double", "c0", "c0")
    c1 = PYB11property("double", "c1", "c1")
    c2pos = PYB11property("double", "c2pos", "c2pos")
    c2neg = PYB11property("double", "c2neg", "c2neg")
    E0 = PYB11property("double", "E0", "E0")
    c0 = PYB11property("double", "c0", "c0")

    Cv = PYB11property("double", "Cv")
    atomicWeight = PYB11property("double", "atomicWeight", "atomicWeight")
    externalPressure = PYB11property("double", "externalPressure", "externalPressure")

#-------------------------------------------------------------------------------
# Inject EOS interface
#-------------------------------------------------------------------------------
PYB11inject(EOSAbstractMethods, OsborneEquationOfState, virtual=True)
