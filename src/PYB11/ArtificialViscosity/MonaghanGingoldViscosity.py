#-------------------------------------------------------------------------------
# MonaghanGingoldViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ArtificialViscosityView import *
from ArtificialViscosityAbstractMethods import *

@PYB11template("Dimension")
class MonaghanGingoldViscosity(ArtificialViscosity):

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
               quadraticInExpansion = ("bool", "false")):
        "MonaghanGingoldViscosity constructor"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"

    #...........................................................................
    # Properties
    linearInExpansion = PYB11property("bool", "linearInExpansion", "linearInExpansion", 
                                      doc="Toggle if the linearviscosity is active for expansion flows")
    quadraticInExpansion = PYB11property("bool", "quadraticInExpansion", "quadraticInExpansion", 
                                         doc="Toggle if the quadratic viscosity is active for expansion flows")

# @PYB11template("Dimension")
# @PYB11template_dict({"QPiType": "typename %(Dimension)s::Scalar"})
# class MonaghanGingoldViscosityView(ArtificialViscosityView):

#     PYB11typedefs = """
#     using Scalar = typename %(Dimension)s::Scalar;
#     using Vector = typename %(Dimension)s::Vector;
#     using Tensor = typename %(Dimension)s::Tensor;
#     using SymTensor = typename %(Dimension)s::SymTensor;
#     using ThirdRankTensor = typename %(Dimension)s::ThirdRankTensor;
#     using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
#     using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
#     using ReturnType = %(QPiType)s;
# """

#     #...........................................................................
#     # Constructors
#     def pyinit(self,
#                Clinear = "const Scalar",
#                Cquadratic = "const Scalar",
#                linearInExpansion = ("bool", "false"),
#                quadraticInExpansion = ("bool", "false")):
#         "MonaghanGingoldViscosityView constructor"    
# #-------------------------------------------------------------------------------
# # Inject abstract interface
# #-------------------------------------------------------------------------------
# PYB11inject(ArtificialViscosityAbstractMethods, MonaghanGingoldViscosityView, virtual=True, pure_virtual=False)
