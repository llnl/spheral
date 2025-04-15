#-------------------------------------------------------------------------------
# PredictorCorrectorIntegrator
#-------------------------------------------------------------------------------
from PYB11Generator import *
from IntegratorAbstractMethods import *
from Integrator import *

@PYB11template("Dimension")
@PYB11cppname("PredictorCorrector")
class PredictorCorrectorIntegrator(Integrator):
    "Second-order in time predictor-corrector integration scheme"

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
PYB11inject(IntegratorAbstractMethods, PredictorCorrectorIntegrator, pure_virtual=False, virtual=True)
