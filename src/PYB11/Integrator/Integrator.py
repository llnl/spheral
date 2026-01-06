#-------------------------------------------------------------------------------
# Integrator base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from IntegratorAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralIntegrator")
class Integrator:
    "Base class for all Spheral time integration algorithms"

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ThirdRankTensor = typename %(Dimension)s::ThirdRankTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               dataBase = "DataBase<%(Dimension)s>&",
               physicsPackages = ("const std::vector<Physics<%(Dimension)s>*>&", "std::vector<Physics<%(Dimension)s>*>()")):
        "Construct an integrator with a DataBase and optional physics packages"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11pycppname("step")
    def step1(self, maxTime="Scalar"):
        "Take a step"
        return "bool"

    @PYB11virtual
    @PYB11const
    def selectDt(self,
                 dtMin = "const Scalar",
                 dtMax = "const Scalar",
                 state = "const State<%(Dimension)s>&",
                 derivs = "const StateDerivatives<%(Dimension)s>&"):
        "Provide a method of looping over the physics packages and picking a time step."
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def preStepInitialize(self,
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        """Perform generic initializations at the beginning of a timestep.
To be called once per advance cycle."""
        return "void"

    @PYB11virtual
    @PYB11const
    def initializeDerivatives(self,
                              t = "const double",
                              dt = "const double",
                              state = "State<%(Dimension)s>&",
                              derivs = "StateDerivatives<%(Dimension)s>&"):
        """Prepare all physics packages for calls to evaluateDerivatives.
To be called before any call to Physics::evaluateDerivatives, therefore potentially
several times during a time step."""
        return "void"

    @PYB11virtual
    @PYB11const
    def postStepFinalize(self,
                         t = "const double",
                         dt = "const double",
                         state = "State<%(Dimension)s>&",
                         derivs = "StateDerivatives<%(Dimension)s>&"):
        "Finalize at the end a timestep, therefore called once at the end of a timestep."
        return "void"

    #...........................................................................
    # Methods
    @PYB11const
    def evaluateDerivatives(self,
                            t = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Iterate over all physics packages and call evaluateDerivatives."
        return "void"

    @PYB11const
    def finalizeDerivatives(self,
                            t = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Iterate over all physics packages and call finalizeDerivatives."
        return "void"

    @PYB11const
    def postStateUpdate(self,
                        t = "const Scalar",
                        dt = "const Scalar",
                        dataBase = "const DataBase<%(Dimension)s>&",
                        state = "State<%(Dimension)s>&",
                        derivs = "StateDerivatives<%(Dimension)s>&"):
        "Iterate over all physics packages and call postStateUpdate"
        return "void"

    def appendPhysicsPackage(self, package="Physics<%(Dimension)s>&"):
        "Add a Physics package."
        return "void"

    def resetPhysicsPackages(self, packages="std::vector<Physics<%(Dimension)s>*>&"):
        "Reset the list of Physics packages."
        return "void"

    @PYB11const
    def havePhysicsPackage(self, package="const Physics<%(Dimension)s>&"):
        "Test if the given Physics package is listed in the integrator."
        return "bool"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def physicsPackages(self):
        "The set of Physics packages registered with the Integrator"
        return "const std::vector<Physics<%(Dimension)s>*>&"

    @PYB11const
    def uniqueBoundaryConditions(self):
        "Get the unique set of boundary conditions across all physics packages."
        return "std::vector<Boundary<%(Dimension)s>*>"

    @PYB11const
    def setGhostNodes(self):
        "Set the ghost nodes for all node lists according to the boundary conditions."
        return "void"

    @PYB11const
    def applyGhostBoundaries(self,
                             state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Set the ghost node values on the Fields of the nodes lists in the data base."
        return "void"

    @PYB11const
    def finalizeGhostBoundaries(self):
        "Finalize the ghost node boundary conditions."
        return "void"

    @PYB11const
    def setViolationNodes(self):
        "Find the nodes in violation of the boundary conditions."
        return "void"

    @PYB11const
    def enforceBoundaries(self,
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        """Reset any internal nodes in violation of boundary conditions to be brought 
into compliance."""
        return "void"

    @PYB11const
    def copyGhostState(self,
                       state0 = "const State<%(Dimension)s>&",
                       state1 = "State<%(Dimension)s>&"):
        "Copy the ghost positions and H's from one state to another."
        return "void"

    #...........................................................................
    # Properties
    currentTime = PYB11property("Scalar", "currentTime", "currentTime", doc="Simulation time")
    currentCycle = PYB11property("int", "currentCycle", "currentCycle", doc="Simulation cycle")
    dtMin = PYB11property("Scalar", "dtMin", "dtMin", doc="Minimum allowed time step")
    dtMax = PYB11property("Scalar", "dtMax", "dtMax", doc="Maximum allowed time step")
    lastDt = PYB11property("Scalar", "lastDt", "lastDt", doc="Last timestep used")
    dtGrowth = PYB11property("Scalar", "dtGrowth", "dtGrowth", doc="Maximum allowed fractional time step growth")
    dtCheckFrac = PYB11property("Scalar", "dtCheckFrac", "dtCheckFrac", doc="The fraction of the timestep we consider when checking for stable behavior")
    dataBase = PYB11property("const DataBase<%(Dimension)s>&", "dataBase", returnpolicy="reference_internal", doc="The DataBase of NodeLists")
    #physicsPackages = PYB11property("const std::vector<Physics<%(Dimension)s>*>&", returnpolicy="reference_internal", doc="The set of physics packages")
    rigorousBoundaries = PYB11property("bool", "rigorousBoundaries", "rigorousBoundaries", doc="Toggle if ghost nodes should be recomputed every derivative estimate")
    updateBoundaryFrequency = PYB11property("int", "updateBoundaryFrequency", "updateBoundaryFrequency", doc="Optionally update the boundary ghost nodes only on this frequency of cycles")
    verbose = PYB11property("bool", "verbose", "verbose", doc="Verbose time step information every step")
    verboseStep = PYB11property("int", "verboseStep", "verboseStep", doc="Select frequency of verbose time step output")
    allowDtCheck = PYB11property("bool", "allowDtCheck", "allowDtCheck", doc="Should the integrator check interim timestep votes and abort steps?")
    domainDecompositionIndependent = PYB11property("bool", "domainDecompositionIndependent", "domainDecompositionIndependent", doc="Order operations to be bit perfect reproducible regardless of domain decomposition")
    cullGhostNodes = PYB11property("bool", "cullGhostNodes", "cullGhostNodes", doc="Cull ghost nodes to just active set")
    dtMultiplier = PYB11property("Scalar", "dtMultiplier", doc="Multiplier for selected time step")

#-------------------------------------------------------------------------------
# Inject other interfaces
#-------------------------------------------------------------------------------
PYB11inject(IntegratorAbstractMethods, Integrator, virtual=False, pure_virtual=True)
PYB11inject(RestartMethods, Integrator)
