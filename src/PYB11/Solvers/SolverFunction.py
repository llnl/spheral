#-------------------------------------------------------------------------------
# SolverFunction
#-------------------------------------------------------------------------------
from PYB11Generator import *

#@PYB11pyname("SolverFunction")
class SolverFunction:
    "Provides base class for functions that compute resdiuals for systems of equation to be solved"

    #...........................................................................
    # Constructors
    def pyinit(self,
               numUnknowns = ("size_t", "0u")):
        return

    #...........................................................................
    # Methods
    @PYB11pure_virtual
    @PYB11const
    def __call__(self,
                 residuals = "std::vector<double>&",
                 x = "const std::vector<double>&"):
        """Override this method to provide the residuals for a given estimate of the system unkonwns.
x and residuals should be of the same length."""
        return "void"

    #...........................................................................
    # Properties
    numUnknowns = PYB11property("size_t", "numUnknowns", "numUnknowns")
