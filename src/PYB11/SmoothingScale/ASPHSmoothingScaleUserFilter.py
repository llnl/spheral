#-------------------------------------------------------------------------------
# ASPHSmoothingScaleUserFilter
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
@PYB11holder("std::shared_ptr")
class ASPHSmoothingScaleUserFilter:

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "ASPHSmoothingScaleUserFilter constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def startFinalize(self,
                      time = "const Scalar", 
                      dt = "const Scalar",
                      dataBase = "DataBase<%(Dimension)s>&", 
                      state = "State<%(Dimension)s>&",
                      derivs = "StateDerivatives<%(Dimension)s>&"):
        "Called at the beginning of ASPHSmoothingScale::finalize"
        return "void"

    @PYB11virtual
    def __call__(self,
                 nodeListi = "size_t",
                 i = "size_t",
                 H0 = "const SymTensor&",
                 H1 = "const SymTensor&"):
        "Called for each point with the old (H0) and new (H1) votes for H(nodeList, i).  Returns the new H to use."
        return "SymTensor"
