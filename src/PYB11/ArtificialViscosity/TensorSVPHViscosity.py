#-------------------------------------------------------------------------------
# TensorSVPHViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ArtificialViscosity import *
from ArtificialViscosityAbstractMethods import *

@PYB11template("Dimension")
@PYB11template_dict({"QPiType": "typename %(Dimension)s::Tensor"})
class TensorSVPHViscosity(ArtificialViscosity):

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
               kernel = "const TableKernel<%(Dimension)s>&",
               fslice = "const Scalar"):
        "TensorSVPHViscosity constructor"

    #...........................................................................
    # Methods
    @PYB11virtual
    def initialize(self,
                   time = "const Scalar", 
                   dt = "const Scalar",
                   dataBase = "const DataBase<%(Dimension)s>&", 
                   state = "State<%(Dimension)s>&",
                   derivs = "StateDerivatives<%(Dimension)s>&"):
        "Some packages might want a hook to do some initializations before the evaluateDerivatives() method is called."
        return "bool"

    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"

    #...........................................................................
    # Properties
    fslice = PYB11property("Scalar", "fslice", "fslice")
    DvDx = PYB11property("const std::vector<Tensor>&", "DvDx", returnpolicy="reference_internal")
    shearCorrection = PYB11property("const std::vector<Scalar>&", "shearCorrection", returnpolicy="reference_internal")
    Qface = PYB11property("const std::vector<Tensor>&", "Qface", returnpolicy="reference_internal")
    
#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(ArtificialViscosityAbstractMethods, TensorSVPHViscosity, virtual=True, pure_virtual=False)
