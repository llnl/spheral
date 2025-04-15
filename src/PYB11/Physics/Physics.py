#-------------------------------------------------------------------------------
# Physics base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from PhysicsAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralPhysics")
class Physics:

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
        "Physics default constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def dtImplicit(self,
                   dataBase = "const DataBase<%(Dimension)s>&",
                   state = "const State<%(Dimension)s>&",
                   derivs = "const StateDerivatives<%(Dimension)s>&",
                   currentTime = "const Scalar"):
        "Optionally compute a timestep for implicit time advancement"
        return "TimeStepType"

    @PYB11virtual
    @PYB11const
    def maxResidual(self,
                    dataBase = "const DataBase<%(Dimension)s>&",
                    state1 = "const State<%(Dimension)s>&",
                    state0 = "const State<%(Dimension)s>&",
                    tol = "const Scalar"):
        "Compute the maximum residual difference between the States"
        return "ResidualType"

    @PYB11virtual
    def applyGhostBoundaries(self,
                             state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    def enforceBoundaries(self,
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Enforce boundary conditions for the physics specific fields."
        return "void"

    @PYB11virtual
    def initializeProblemStartup(self,
                                 dataBase = "DataBase<%(Dimension)s>&"):
        """An optional hook to initialize once when the problem is starting up.
Typically this is used to size arrays once all the materials and NodeLists have
been created.  It is assumed after this method has been called it is safe to
call Physics::registerState for instance to create full populated State objects."""

        return "void"

    @PYB11virtual
    def initializeProblemStartupDependencies(self,
                                             dataBase = "DataBase<%(Dimension)s>&",
                                             state = "State<%(Dimension)s>&",
                                             derivs = "StateDerivatives<%(Dimension)s>&"):
        """A second optional method to be called on startup, after Physics::initializeProblemStartup has
been called.
One use for this hook is to fill in dependendent state using the State object, such as
temperature or pressure."""
        return "void"

    @PYB11virtual
    def preStepInitialize(self,
                          dataBase = "const DataBase<%(Dimension)s>&", 
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Optional hook to be called at the beginning of a time step."
        return "void"

    @PYB11virtual
    def initialize(self,
                   time = "const Scalar", 
                   dt = "const Scalar",
                   dataBase = "const DataBase<%(Dimension)s>&", 
                   state = "State<%(Dimension)s>&",
                   derivs = "StateDerivatives<%(Dimension)s>&"):
        "Some packages might want a hook to do some initializations before the evaluateDerivatives() method is called.  Returns a bool indicating if boundary conditions should be applied."
        return "bool"

    @PYB11virtual
    def finalize(self,
                 time = "const Scalar", 
                 dt = "const Scalar",
                 dataBase = "DataBase<%(Dimension)s>&", 
                 state = "State<%(Dimension)s>&",
                 derivs = "StateDerivatives<%(Dimension)s>&"):
        "Similarly packages might want a hook to do some post-step finalizations.  Really we should rename this post-step finalize."
        return "void"

    @PYB11virtual
    @PYB11const
    def finalizeDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Provide a hook to be called after all physics packages have had their evaluateDerivatives method called, but before anyone does anything with those derivatives."
        return "void"

    @PYB11virtual
    def postStateUpdate(self,
                        time = "const Scalar", 
                        dt = "const Scalar",
                        dataBase = "const DataBase<%(Dimension)s>&", 
                        state = "State<%(Dimension)s>&",
                        derivs = "StateDerivatives<%(Dimension)s>&"):
        """Provide a hook to be called after the state has been updated and boundary conditions have been enforced.
Returns a bool indicating if boundary conditions should be applied."""
        return "bool"

    @PYB11virtual
    @PYB11const
    def requireConnectivity(self):
        "Some physics does not require the connectivity be constructed."
        return "bool"

    @PYB11virtual
    @PYB11const
    def requireGhostConnectivity(self):
        "Some physics algorithms require ghost connectivity to be constructed."
        return "bool"

    @PYB11virtual
    @PYB11const
    def requireOverlapConnectivity(self):
        "Some physics algorithms require overlap connectivity to be constructed."
        return "bool"

    @PYB11virtual
    @PYB11const
    def requireIntersectionConnectivity(self):
        "Some physics algorithms require intersection connectivity to be constructed."
        return "bool"

    @PYB11virtual
    @PYB11const
    def requireVoronoiCells(self):
        "Some physics algorithms require the Voronoi cells per point be computed."
        return "bool"

    @PYB11virtual
    @PYB11const
    def requireReproducingKernels(self):
        "Some physics algorithms require reproducing kernels."
        return "std::set<RKOrder>"

    @PYB11virtual
    @PYB11const
    def requireReproducingKernelHessian(self):
        "Some physics algorithms require reproducing kernel spatial second derivatives."
        return "bool"

    @PYB11virtual
    @PYB11const
    def updateReproducingKernelsInFinalize(self):
        "Does this package need an update of reproducing kernels during finalize?"
        return "bool"
    
    @PYB11virtual
    @PYB11const
    def extraEnergy(self):
        "Many physics packages will have their own representations of energy in the system (gravitational potential energy, radiative losses, etc.)"
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def extraMomentum(self):
        "Many physics packages will also have their own representations of momentum in the system (electromagnetic momentum flux density, etc.)"
        return "Vector"

    @PYB11virtual
    def registerAdditionalVisualizationState(self,
                                             dataBase = "DataBase<%(Dimension)s>&",
                                             state = "State<%(Dimension)s>&"):
        "Register any additional state for visualization."
        return "void"

    #...........................................................................
    # Methods
    def appendBoundary(self, boundary="Boundary<%(Dimension)s>&"):
        "Add a Boundary condition to the end of boundary list"
        return "void"
    
    def prependBoundary(self, boundary="Boundary<%(Dimension)s>&"):
        "Insert a Boundary condition at the beginning of the boundary list"
        return "void"

    def clearBoundaries(self):
        "Remove all boundary conditions"
        return "void"

    @PYB11const
    def haveBoundary(self, boundary="const Boundary<%(Dimension)s>&"):
        "Test if the given Boundary condition is registered."
        return "bool"

    # @PYB11returnpolicy("reference_internal")
    # @PYB11const
    # def boundaryConditions(self):
    #     "Access the list of boundary conditions."
    #     return "const std::vector<Boundary<%(Dimension)s>*>&"

    def appendSubPackage(self, package="Physics<%(Dimension)s>&"):
        "Add a package to be run after this one"
        return "void"
    
    def prependSubPackage(self, package="Physics<%(Dimension)s>&"):
        "Add a package to run before this one"
        return "void"

    #...........................................................................
    # Properties
    #"std::vector<Boundary<%(Dimension)s>*>", 
    boundaryConditions = PYB11property(doc="The set of boundary conditions")
    postSubPackages = PYB11property(doc="Packages that should be run after this one")
    preSubPackages = PYB11property(doc="Packages that should be run before this one")

#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, Physics, pure_virtual=True)
