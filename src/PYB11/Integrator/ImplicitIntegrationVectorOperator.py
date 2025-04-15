#-------------------------------------------------------------------------------
# ImplicitIntegrationVectorOperator
#-------------------------------------------------------------------------------
from PYB11Generator import *
import SolverFunction

@PYB11template("Dimension")
class ImplicitIntegrationVectorOperator(SolverFunction.SolverFunction):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               t = "const Scalar",
               dt = "const Scalar",
               beta = "const Scalar",
               state0 = "const State<%(Dimension)s>&",
               derivs0 = "const StateDerivatives<%(Dimension)s>&",
               state = "State<%(Dimension)s>&",
               derivs = "StateDerivatives<%(Dimension)s>&",
               integrator = "Integrator<%(Dimension)s>&"):
        return

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def __call__(self,
                 residuals = "std::vector<double>&",
                 x = "const std::vector<double>&"):
        return "void"

    #...........................................................................
    # Properties
    inputState0 = PYB11property("const std::vector<double>&", "inputState0")
    derivs0 = PYB11property("const std::vector<double>&", "derivs0")
    state = PYB11property("State<%(Dimension)s>&", "state", "state")
    derivs = PYB11property("StateDerivatives<%(Dimension)s>&", "derivs", "derivs")
    integrator = PYB11property("Integrator<%(Dimension)s>&", "integrator", "integrator")
    t = PYB11property("Scalar", "t", "t")
    dt = PYB11property("Scalar", "dt", "dt")
    beta = PYB11property("Scalar", "beta", "beta")
