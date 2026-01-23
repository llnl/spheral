#-------------------------------------------------------------------------------
# GenericBodyForce base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralPhysics")
class GenericBodyForce(Physics):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ThirdRankTensor = typename %(Dimension)s::ThirdRankTensor;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "GenericBodyForce constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&", 
                      state = "State<%(Dimension)s>&"):
        "Default state registration for an acceleration source"
        return "void"

    @PYB11virtual
    def registerDerivatives(self,
                            dataBase = "DataBase<%(Dimension)s>&", 
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Default state derivative registration for an acceleration source"
        return "void"

    @PYB11pure_virtual
    @PYB11const
    def dt(dataBase = "const DataBase<%(Dimension)s>&", 
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"
    #...........................................................................
    # Attributes
    DxDt = PYB11property("const FieldList<%(Dimension)s, Vector>&", "DxDt", doc="Time derivative for position")
    DvDt = PYB11property("const FieldList<%(Dimension)s, Vector>&", "DvDt", doc="Time derivative for velocity")
