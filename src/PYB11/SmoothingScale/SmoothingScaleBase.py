#-------------------------------------------------------------------------------
# SmoothingScaleBase abstract class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *

@PYB11template("Dimension")
class SmoothingScaleBase(Physics):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ThirdRankTensor = typename %(Dimension)s::ThirdRankTensor;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
    using HidealFilterType = typename ASPHSmoothingScale<%(Dimension)s>::HidealFilterType;
    using RadialFunctorType = typename ASPHSmoothingScale<%(Dimension)s>::RadialFunctorType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               HUpdate = "HEvolutionType",
               fixShape = "const bool",
               radialOnly = "const bool"):
        "SmoothingScaleBase constructor"

    #...........................................................................
    # Pure virtual methods
    @PYB11pure_virtual
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Increment the derivatives."
        return "void"

    @PYB11pure_virtual
    @PYB11const
    def label(self):
        "It's useful to have labels for Physics packages.  We'll require this to have the same signature as the restart label."
        return "std::string"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def initializeProblemStartup(self,
                                 dataBase = "DataBase<%(Dimension)s>&"):
        """An optional hook to initialize once when the problem is starting up.
Typically this is used to size arrays once all the materials and NodeLists have
been created.  It is assumed after this method has been called it is safe to
call Physics::registerState for instance to create full populated State objects."""
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
    def preStepInitialize(self,
                          dataBase = "const DataBase<%(Dimension)s>&", 
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Optional hook to be called at the beginning of a time step."
        return "void"

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
    @PYB11const
    def dumpState(self, file="FileIO&", pathName="const std::string&"):
        "Serialize under the given path in a FileIO object"
        return "void"

    @PYB11virtual
    def restoreState(self, file="const FileIO&", pathName="const std::string&"):
        "Restore state from the given path in a FileIO object"
        return "void"

    #...........................................................................
    # Methods
    @PYB11const
    def hmax(self,
             Vi = "const Scalar",
             nPerh = "const Scalar"):
        "Given the volume and target nperh, compute an effective target hmax"
        return "Scalar"

    #...........................................................................
    # Attributes
    HEvolution = PYB11property("HEvolutionType", "HEvolution", "HEvolution", doc="The H evolution choice")
    Hideal = PYB11property("const FieldList<%(Dimension)s, SymTensor>&", "Hideal", doc="The ideal H storage FieldList")
    DHDt = PYB11property("const FieldList<%(Dimension)s, SymTensor>&", "DHDt", doc="The H time derivative storage FieldList")
    radius0 = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "radius0", doc="The radius of a point at the beginning of a timestep (if using radialOnly=True)")
    HidealFilter = PYB11property("std::shared_ptr<HidealFilterType>", "HidealFilter", "HidealFilter", doc="Optional functor to manipulate the Hideal calculation")
    RadialFunctor = PYB11property("std::shared_ptr<RadialFunctorType>", "RadialFunctor", "RadialFunctor", doc="Optional functor to manipulate the radial normal and radius are computed when using radialOnly")
    damping = PYB11property("double", "damping", "damping", doc="exponent on the change in smoothing parameter, e.g. 0.5 / dim, for use in iterateIdealH")
    fixShape = PYB11property("bool", "fixShape", "fixShape", doc="Force the H tensor shape to be fixed -- only adjust volume")
    radialOnly = PYB11property("bool", "radialOnly", "radialOnly", doc="Force the H tensor to evolve solely in the radial direction")
