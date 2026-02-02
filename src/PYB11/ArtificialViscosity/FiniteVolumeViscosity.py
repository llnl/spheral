#-------------------------------------------------------------------------------
# FiniteVolumeViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ArtificialViscosity import *
from ArtificialViscosityAbstractMethods import *

@PYB11template("Dimension")
@PYB11template_dict({"QPiType": "typename %(Dimension)s::Scalar"})
class FiniteVolumeViscosity(ArtificialViscosity):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ThirdRankTensor = typename %(Dimension)s::ThirdRankTensor;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
    using ReturnType = %(QPiType)s;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               Clinear = "const Scalar",
               Cquadratic = "const Scalar",
               kernel = "const TableKernel<%(Dimension)s>&"):
        "FiniteVolumeViscosity constructor"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def requireVelocityGradient(self):
        "We need the velocity gradient and set this to true"
        return "bool"

    @PYB11virtual
    def updateVelocityGradient(self,
                               dataBase = "const DataBase<%(Dimension)s>&",
                               state = "const State<%(Dimension)s>&",
                               derivs = "const StateDerivatives<%(Dimension)s>&"):
        "Update the locally stored velocity gradient"
        return "void"

    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"
    
#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(ArtificialViscosityAbstractMethods, FiniteVolumeViscosity, virtual=True, pure_virtual=False)
