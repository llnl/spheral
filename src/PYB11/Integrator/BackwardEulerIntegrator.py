#-------------------------------------------------------------------------------
# BackwardEulerIntegrator
#-------------------------------------------------------------------------------
from PYB11Generator import *
from IntegratorAbstractMethods import *
from ImplicitIntegrator import *

@PYB11template("Dimension")
@PYB11cppname("BackwardEuler")
class BackwardEulerIntegrator(ImplicitIntegrator):
    "First-order in time implicit (backward Euler) integration scheme"

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
               packages = ("vector<Physics<%(Dimension)s>*>", "vector<Physics<%(Dimension)s>*>()"),
               beta = ("Scalar", "1.0"),
               ftol = ("Scalar", "1.0e-8"),
               tol = ("Scalar", "1.0e-8"),
               maxIterations = ("size_t", "100u")):
        """Construct a backward Euler itegrator.
Note: beta is the blending of the (n+1) and (n) time derivatives, so adjusting this parameter
makes this not just backward Euler:
    beta = 1.0 : backward Euler
    beta = 0.0 : forward Euler
    beta = 0.5 : Crank-Nicholson"""
        return

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def step(self,
             maxTime = "Scalar",
             state = "State<%(Dimension)s>&",
             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Take a step"
        return "bool"

    #...........................................................................
    # Properties
    beta = PYB11property("Scalar", "beta", "beta", doc="The blend of (n+1) and (n) derivative states for advancement")
    maxIterations = PYB11property("size_t", "maxIterations", "maxIterations", doc="The maximum allowed iterations to try for advancing a step")
    ftol = PYB11property("Scalar", "ftol", "ftol")
    tM2 = PYB11property("Scalar", "tM2")
    tM1 = PYB11property("Scalar", "tM2")
    solutionM2 = PYB11property("const std::vector<double>&", "solutionM2")
    solutionM1 = PYB11property("const std::vector<double>&", "solutionM2")

#-------------------------------------------------------------------------------
# Inject other interfaces
#-------------------------------------------------------------------------------
PYB11inject(IntegratorAbstractMethods, BackwardEulerIntegrator, pure_virtual=False, virtual=True)
