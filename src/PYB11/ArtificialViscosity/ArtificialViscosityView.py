#-------------------------------------------------------------------------------
# ArtificialViscosityView
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ArtificialViscosity import *
from ArtificialViscosityAbstractMethods import *

@PYB11template("Dimension", "QPiType")
@PYB11module("SpheralArtificialViscosity")
class ArtificialViscosityView():

    PYB11typedefs = """
  using Scalar = typename %(Dimension)s::Scalar;
  using Vector = typename %(Dimension)s::Vector;
  using Tensor = typename %(Dimension)s::Tensor;
  using SymTensor = typename %(Dimension)s::SymTensor;
  using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
  using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
  using ReturnType = %(QPiType)s;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               Clinear = "const Scalar",
               Cquadratic = "const Scalar"):
        "ArtificialViscosityView constructor"

    #...........................................................................
    # Methods
    @PYB11const
    def QPiTypeIndex(self):
        "Require ArtificialViscosities to specify the type_index of the descendant QPiType"
        return "std::type_index"
    
#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(ArtificialViscosityAbstractMethods, ArtificialViscosityView, pure_virtual=True)
