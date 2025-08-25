#-------------------------------------------------------------------------------
# TensorMonaghanGingoldViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ArtificialViscosityView import *
from ArtificialViscosityAbstractMethods import *

@PYB11template("Dimension")
@PYB11template_dict({"QPiType": "typename %(Dimension)s::Tensor"})
class TensorMonaghanGingoldViscosity(ArtificialViscosityView):
    """A modified form of the Monaghan & Gingold viscosity, extended to tensor formalism.
This method is described in
Owen, J Michael (2004), 'A tensor artficial visocity for SPH', Journal of Computational Physics, 201(2), 601-629
"""

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
        "TensorMonaghanGingoldViscosity constructor"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def requireVelocityGradient(self):
        "We need the velocity gradient and set this to true"
        return "bool"

    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"
    
#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(ArtificialViscosityAbstractMethods, TensorMonaghanGingoldViscosity, virtual=True, pure_virtual=False)
