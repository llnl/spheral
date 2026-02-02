#-------------------------------------------------------------------------------
# PorosityModel
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralPorosity")
class PorosityModel(Physics):
    """PorosityModel
Base class for PorosityModels for common functionality.
"""

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using ScalarField = Field<%(Dimension)s, Scalar>;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               nodeList = "const SolidNodeList<%(Dimension)s>&",
               phi0 = "const double",
               cS0 = "const double",
               c0 = "const double",
               rhoS0 = "const double",
               jutziStateUpdate = "bool"):
        """Constructor parameters:
        nodeList:         The SolidNodeList we're going apply to
        phi0:             Initial porosity (single value)
        cS0:              Reference sound speed at full density
        c0:               Reference sound speed at initial porosity
        rhoS0:            Reference solid density
        jutziStateUpdate: update state as describd in Jutzi 2008"""

    def pyinit1(self,
                nodeList = "const SolidNodeList<%(Dimension)s>&",
                phi0 = "const Field<%(Dimension)s, %(Dimension)s::Scalar>&",
                cS0 = "const double",
                c0 = "const Field<%(Dimension)s, %(Dimension)s::Scalar>&",
                rhoS0 = "const double",
                jutziStateUpdate = "bool"):
        """Constructor parameters:
        nodeList:         The SolidNodeList we're going apply to
        phi0:             Initial porosity (field of values)
        cS0:              Reference sound speed at full density
        c0:               Reference sound speed at initial porosity (field)
        rhoS0:            Reference solid density
        jutziStateUpdate: update state as describd in Jutzi 2008"""

    #...........................................................................
    # Virtual methods
    @PYB11const
    @PYB11virtual
    def dt(dataBase = "const DataBase<%(Dimension)s>&", 
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

    @PYB11virtual
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state you want carried around (and potentially evolved), as well as the policies for such evolution."
        return "void"

    @PYB11virtual
    def registerDerivatives(self,
                            dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Register the derivatives/change fields for updating state."
        return "void"


    @PYB11virtual
    def initializeProblemStartup(self, 
                                 dataBase = "DataBase<%(Dimension)s>&"):
        "Do any required one-time initializations on problem start up."
        return "void"

    #...........................................................................
    # Methods
    @PYB11const
    def phi(self):
        "Compute the current porosity from the distention"
        return "ScalarField"

    #...........................................................................
    # Properties
    jutziStateUpdate = PYB11property(doc="Switch to update state as described in Jutzi 2008")
    rhoS0 = PYB11property(doc="Reference solid density")
    cS0 = PYB11property(doc="Reference sound speed at full density")
    KS0 = PYB11property()
    fdt = PYB11property("double", getter="fdt", setter="fdt", doc="The timestep fractional multiplier (0 => no timestep control on alpha)")
    maxAbsDalphaDt = PYB11property(doc="maximum of the last abs(DalphaDt) calculated")
    nodeList = PYB11property("const SolidNodeList<%(Dimension)s>&", returnpolicy="reference_internal")
    alpha0 = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    alpha = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    DalphaDt = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    solidMassDensity = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    c0 = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    fDS = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    fDSnew = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
