#-------------------------------------------------------------------------------
# CheapSynchronousRK2Integrator
#-------------------------------------------------------------------------------
from PYB11Generator import *
from IntegratorAbstractMethods import *
from Integrator import *

@PYB11template("Dimension")
@PYB11cppname("CheapSynchronousRK2")
class CheapSynchronousRK2Integrator(Integrator):
    """Second-order in time explicit Runge-Kutta time integration scheme
This method cheats and reuses the previous step mid-timestep time derivatives
to predict the state at the middle of the current step, rather than evaluate
the derivatives at the beginning of this step.  This allows us to only do a single
call to evaluate derivatives per timestep, maintaining formal second-order but 
sacrificing some accuracy vs. the the true RK2 algorithm."""

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

#-------------------------------------------------------------------------------
# Inject other interfaces
#-------------------------------------------------------------------------------
PYB11inject(IntegratorAbstractMethods, CheapSynchronousRK2Integrator, pure_virtual=False, virtual=True)
