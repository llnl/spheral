#-------------------------------------------------------------------------------
# PolyGravity base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from RestartMethods import *
from GenericBodyForce import *

@PYB11template("Dimension")
class PolyGravity(GenericBodyForce):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ThirdRankTensor = typename %(Dimension)s::ThirdRankTensor;
    using Polytope = typename %(Dimension)s::FacetedVolume;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               poly = "const Polytope&",
               G = "const double",
               mass = "const double",
               ftimestep = "const double",
               timeStepChoice = ("GravityTimeStepType", "GravityTimeStepType::DynamicalTime")):
        "PolyGravity constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state you want carried around (and potentially evolved), as well as the policies for such evolution."
        return "void"

    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Increment the derivatives."
        return "void"

    @PYB11virtual
    @PYB11const
    def dt(dataBase = "const DataBase<%(Dimension)s>&", 
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

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
    @PYB11const
    def requireConnectivity(self):
        "Some physics does not require the connectivity be constructed."
        return "bool"

    @PYB11virtual
    @PYB11const
    def extraEnergy(self):
        "Many physics packages will have their own representations of energy in the system (gravitational potential energy, radiative losses, etc.)"
        return "Scalar"

    #...........................................................................
    # Properties
    potential = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "potential", returnpolicy="reference_internal", doc="The last computed potential")
    poly = PYB11property("const Dim<3>::FacetedVolume&", "poly", returnpolicy="reference_internal", doc="The 3D polyhedral surface")
    G = PYB11property("double", "G", doc="The gravitational constant")
    mass = PYB11property("double", "mass", doc="The mass of the polyhedron")
    ftimestep = PYB11property("double", "ftimestep", "ftimestep", doc="The current time step scaling factor.")
    timeStepChoice = PYB11property("GravityTimeStepType", "timeStepChoice", "timeStepChoice", doc="The algorithmic choice for setting the time step.")
    dynamicalTime = PYB11property("double", "dynamicalTime", doc="The dynamical time (1/sqrt(G*rho))")
    solver = PYB11property("const ApproximatePolyhedralGravityModel&", "solver", returnpolicy="reference_internal", doc="The 3D gravity solver")

#-------------------------------------------------------------------------------
# Inject restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, PolyGravity)
