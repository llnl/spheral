#-------------------------------------------------------------------------------
# LimitedMonaghanGingoldViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from MonaghanGingoldViscosity import *
from ArtificialViscosityAbstractMethods import *

@PYB11template("Dimension")
class LimitedMonaghanGingoldViscosity(MonaghanGingoldViscosity):

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
               Clinear = "const Scalar",
               Cquadratic = "const Scalar",
               kernel = "const TableKernel<%(Dimension)s>&",
               linearInExpansion = ("bool", "false"),
               quadraticInExpansion = ("bool", "false"),
               etaCritFrac = ("double", "1.0"),
               etaFoldFrac = ("double", "0.2")):
        "LimitedMonaghanGingoldViscosity constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def requireVelocityGradient(self):
        "We need the velocity gradient and set this to true"
        return "bool"

    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"

    #...........................................................................
    # Properties
    etaCritFrac = PYB11property("double", "etaCritFrac", "etaCritFrac",
                                doc="Critical eta, below which viscosity becomes active regardless of limiter")
    etaFoldFrac = PYB11property("double", "etaFoldFrac", "etaFoldFrac",
                                doc="Smoothness in eta phasing in viscosity around etaCritFrac")

@PYB11template("Dimension")
@PYB11template_dict({"QPiType": "typename %(Dimension)s::Scalar"})
class LimitedMonaghanGingoldViscosityView(MonaghanGingoldViscosityView):

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
               linearInExpansion = ("bool", "false"),
               quadraticInExpansion = ("bool", "false"),
               etaCritFrac = ("double", "1.0"),
               etaFoldFrac = ("double", "0.2")):
        "LimitedMonaghanGingoldViscosityView constructor"

#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(ArtificialViscosityAbstractMethods, LimitedMonaghanGingoldViscosityView, virtual=True, pure_virtual=False)
