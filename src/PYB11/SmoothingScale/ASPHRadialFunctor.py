#-------------------------------------------------------------------------------
# ASPHRadialFunctor
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
@PYB11holder("std::shared_ptr")
class ASPHRadialFunctor:

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "ASPHRadialFunctor constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def radialUnitVector(self,
                         nodeListi = "const size_t",
                         i = "const size_t",
                         posi = "const Vector&"):
        return "Vector"

    @PYB11virtual
    @PYB11const
    def radialCoordinate(self,
                         nodeListi = "const size_t",
                         i = "const size_t",
                         posi = "const Vector&"):
        return "Scalar"
