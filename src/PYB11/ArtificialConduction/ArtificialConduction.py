#-------------------------------------------------------------------------------
# ArtificialConduction
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *

@PYB11template("Dimension")
class ArtificialConduction(Physics):
    "ArtificialConduction -- Artificial smoothing of energy discontinuities"

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
               W = "const TableKernel<%(Dimension)s>&",
               alphaArCond = "const Scalar",
               ACcorrectionOrder = ("const RKOrder", "RKOrder::LinearOrder")):
        "Constructor"
            
        #...........................................................................
        # Virtual methods
        @PYB11virtual
        def initializeProblemStartup(self,
                                     dataBase = "DataBase<%(Dimension)s>&"):
            "Do any required one-time initializations on problem start up."
            return "void"
            
        #...........................................................................
        # Properties
        ACcorrectionOrder = PYB11property("RKOrder", "ACcorrectionOrder", "ACcorrectionOrder",
                                          doc="The correction order")

#-------------------------------------------------------------------------------
# Inject physics interface
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, ArtificialConduction, virtual=True, pure_virtual=False)
