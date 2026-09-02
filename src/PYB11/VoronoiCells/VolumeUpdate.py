#-------------------------------------------------------------------------------
# VolumeUpdate
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11dynamic_attr
class VolumeUpdate(Physics):
    "Base class for volume computation packages"

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
    using VolumeRequirements = typename Physics<%(Dimension)s>::VolumeRequirements;
    using RKRequirements = typename Physics<%(Dimension)s>::RKRequirements;
    using ConnectivityRequirements = typename Physics<%(Dimension)s>::ConnectivityRequirements;
"""

    def pyinit(self,
               volumeType = "const VolumeType",
               W = "const TableKernel<%(Dimension)s>&",
               updateInStep = ("const bool", "true"),
               updateInFinalize = ("const bool", "false")):
        "VolumeUpdate constructor"
        return

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def initializeProblemStartup(self,
                                 dataBase = "DataBase<%(Dimension)s>&"):
        "Allocate volume fields on problem startup."
        return "void"

    @PYB11virtual
    def initializeProblemStartupDependencies(self,
                                             dataBase = "DataBase<%(Dimension)s>&",
                                             state = "State<%(Dimension)s>&",
                                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Compute initial volumes."
        return "void"

    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "No derivatives to evaluate."
        return "void"

    @PYB11virtual
    @PYB11const
    def dt(dataBase = "const DataBase<%(Dimension)s>&",
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "No time step vote."
        return "TimeStepType"

    @PYB11virtual
    def registerState(dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register volume in state."
        return "void"

    @PYB11virtual
    def registerDerivatives(dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "No derivatives to register."
        return "void"

    @PYB11virtual
    def applyGhostBoundaries(state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to volume."
        return "void"

    @PYB11virtual
    def enforceBoundaries(state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Enforce boundary conditions on volume."
        return "void"

    #...........................................................................
    # Properties
    volumeType =     PYB11property("VolumeType", "volumeType", doc="The volume computation method")
    updateInStep =   PYB11property("bool", "updateInStep", doc="Whether to update volumes during explicit step")
    updateInFinalize = PYB11property("bool", "updateInFinalize", doc="Whether to update volumes during implicit finalize")
    volume =         PYB11property("const FieldList<%(Dimension)s, Scalar>&", "volume", returnpolicy="reference_internal")
    volume3d = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "volume3d", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, VolumeUpdate)
