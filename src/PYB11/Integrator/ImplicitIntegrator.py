#-------------------------------------------------------------------------------
# ImplicitIntegrator base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from IntegratorAbstractMethods import *
from Integrator import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralIntegrator")
class ImplicitIntegrator(Integrator):
    "Base class for all Spheral implicit in time integration algorithms"

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
               physicsPackages = ("const std::vector<Physics<%(Dimension)s>*>&", "std::vector<Physics<%(Dimension)s>*>()"),
               tol = ("Scalar", "1.0e-6")):
        "Construct an ImplicitIntegrator with a DataBase, physics packages, and tolerance"

    #...........................................................................
    # Methods
    @PYB11virtual
    def step(self,
             maxTime = "Scalar"):
        "Take an implicit step"
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
    def computeResiduals(self,
                         state1 = "const State<%(Dimension)s>&",
                         state0 = "const State<%(Dimension)s>&",
                         forceVerbose = ("bool", "false")):
        "Compute the maximum residual difference across all physics packages between two States"
        return "Scalar"

    #...........................................................................
    # Properties
    convergenceTolerance = PYB11property("Scalar", "convergenceTolerance", "convergenceTolerance", doc="Tolerance for convergence in residuals during integration")
    maxAllowedDtMultiplier = PYB11property("Scalar", "maxAllowedDtMultiplier", "maxAllowedDtMultiplier")
    maxGoodDtMultiplier = PYB11property("Scalar", "maxGoodDtMultiplier")
    numExplicitSteps = PYB11property("size_t", "numExplicitSteps")
    numImplicitSteps = PYB11property("size_t", "numImplicitSteps")

#-------------------------------------------------------------------------------
# Inject other interfaces
#-------------------------------------------------------------------------------
PYB11inject(IntegratorAbstractMethods, ImplicitIntegrator, virtual=False, pure_virtual=True)
