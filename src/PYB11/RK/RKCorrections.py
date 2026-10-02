#-------------------------------------------------------------------------------
# RKCorrections
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class RKCorrections(Physics):
    "Computes RK correction terms"
    
    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using FacetedVolume = typename %(Dimension)s::FacetedVolume;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
    using VolumeRequirements = typename Physics<%(Dimension)s>::VolumeRequirements;
    using RKRequirements = typename Physics<%(Dimension)s>::RKRequirements;
    using ConnectivityRequirements = typename Physics<%(Dimension)s>::ConnectivityRequirements;
"""

    def pyinit(self,
               orders = "const std::set<RKOrder>",
               dataBase = "const DataBase<%(Dimension)s>&",
               W = "const TableKernel<%(Dimension)s>&",
               needHessian = "const bool",
               updateInStep = ("const bool", "true"),
               updateInFinalize = ("const bool", "false")):
        "Constructor"
        
    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Evaluate derivatives"
        return "void"

    @PYB11virtual 
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state"
        return "void"

    @PYB11virtual
    def registerDerivatives(self,
                            dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Register derivatives"
        return "void"

    @PYB11virtual
    def applyGhostBoundaries(self,
                             state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields"
        return "void"

    @PYB11virtual
    def enforceBoundaries(self,
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Enforce boundary conditions for the physics specific fields"
        return "void"
    
    @PYB11virtual
    def initializeProblemStartup(self,
                                 dataBase = "DataBase<%(Dimension)s>&"):
        "Tasks we do once on problem startup"
        return "void"

    @PYB11virtual
    def initializeProblemStartupDependencies(self,
                                             dataBase = "DataBase<%(Dimension)s>&",
                                             state = "State<%(Dimension)s>&",
                                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Tasks we do once on problem startup"
        return "void"

    @PYB11virtual
    def initialize(self,
                   time = "const Scalar",
                   dt = "const Scalar",
                   dataBase = "const DataBase<%(Dimension)s>&",
                   state = "State<%(Dimension)s>&",
                   derivs = "StateDerivatives<%(Dimension)s>&"):
        "Initialize the Hydro before we start a derivative evaluation."
        return "bool"
                  
    @PYB11virtual
    def finalize(self,
                 time = "const Scalar",
                 dt = "const Scalar",
                 dataBase = "DataBase<%(Dimension)s>&",
                 state = "State<%(Dimension)s>&",
                 derivs = "StateDerivatives<%(Dimension)s>&"):
        "Finalize — recompute corrections at end of step."
        return "void"
                  
    @PYB11const
    @PYB11returnpolicy("reference_internal")
    @PYB11keepalive(0,1)
    def WR(self,
           order = "const RKOrder"):
        "Look up the ReproducingKernel for the given order"
        return "const ReproducingKernel<%(Dimension)s>&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    @PYB11keepalive(0,1)
    def corrections(self,
                    order = "const RKOrder"):
        "Look up the corrections for the given order"
        return "const FieldList<%(Dimension)s, RKCoefficients<%(Dimension)s>>&"

    #...........................................................................
    # Properties
    correctionOrders = PYB11property(doc="The set of spatial orders for the reproducing kernel corrections")
    needHessian = PYB11property(doc="Whether the RK hessian is needed")
    updateInStep = PYB11property(doc="Whether to update corrections during explicit step")
    updateInFinalize = PYB11property(doc="Whether to update corrections during implicit finalize")

    surfaceArea = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "surfaceArea", returnpolicy="reference_internal")
    normal = PYB11property("const FieldList<%(Dimension)s, Vector>&", "normal", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, RKCorrections, virtual=True)
PYB11inject(RestartMethods, RKCorrections)
