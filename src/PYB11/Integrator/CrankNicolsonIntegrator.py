#-------------------------------------------------------------------------------
# CrankNicolsonIntegrator
#-------------------------------------------------------------------------------
from PYB11Generator import *
from IntegratorAbstractMethods import *
from ImplicitIntegrator import *

@PYB11template("Dimension")
@PYB11cppname("CrankNicolson")
class CrankNicolsonIntegrator(ImplicitIntegrator):
    "Second-order in time implicit Crank-Nicolson integration scheme"

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
               beta = ("Scalar", "0.5"),
               tol = ("Scalar", "1.0e-4"),
               maxIterations = ("size_t", "10u")):
        """Construct a Crank-Nicolson itegrator.
Note: beta is the blending of the (n+1) and (n) time solutions in the fixed-point iteration."""
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
    beta = PYB11property("Scalar", "beta", "beta", doc="The blend of (n+1) and (n) state during fixed-point iteration")
    maxIterations = PYB11property("size_t", "maxIterations", "maxIterations", doc="The maximum allowed iterations to try for advancing a step")

#-------------------------------------------------------------------------------
# Inject other interfaces
#-------------------------------------------------------------------------------
PYB11inject(IntegratorAbstractMethods, CrankNicolsonIntegrator, pure_virtual=False, virtual=True)
