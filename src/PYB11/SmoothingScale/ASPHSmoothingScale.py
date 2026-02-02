#-------------------------------------------------------------------------------
# ASPHSmoothingScale
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SmoothingScaleBase import *

@PYB11template("Dimension")
class ASPHSmoothingScale(SmoothingScaleBase):

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
    def pyinit(self,
               HUpdate = "HEvolutionType",
               W = "const TableKernel<%(Dimension)s>&",
               fixShape = ("const bool", "false"),
               radialOnly = ("const bool", "false")):
        "ASPHSmoothingScale constructor"

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
    def finalize(self,
                 time = "const Scalar", 
                 dt = "const Scalar",
                 dataBase = "DataBase<%(Dimension)s>&", 
                 state = "State<%(Dimension)s>&",
                 derivs = "StateDerivatives<%(Dimension)s>&"):
        "Similarly packages might want a hook to do some post-step finalizations.  Really we should rename this post-step finalize."
        return "void"

    @PYB11virtual
    def applyGhostBoundaries(self,
                             state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    @PYB11const
    def requireVoronoiCells(self):
        "Some physics algorithms require the Voronoi cells per point be computed."
        return "bool"

    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"

    #...........................................................................
    # Attributes
    WT = PYB11property("const TableKernel<%(Dimension)s>&", "WT", doc="The interpolation kernel")
    zerothMoment = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "zerothMoment", doc="The zeroth moment storage FieldList")
    secondMoment = PYB11property("const FieldList<%(Dimension)s, SymTensor>&", "secondMoment", doc="The second moment storage FieldList")
    cellSecondMoment = PYB11property("const FieldList<%(Dimension)s, SymTensor>&", "cellSecondMoment", doc="The second moment of the Voronoi cells")
