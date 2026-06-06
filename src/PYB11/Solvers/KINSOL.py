#-------------------------------------------------------------------------------
# KINSOL
#-------------------------------------------------------------------------------
from PYB11Generator import *

class KINSOL:
    "Wrapper around the Sundials KINSOL nonlinear solver"

    #...........................................................................
    # Constructors
    def pyinit(self):
        "KINSOL constructor"
        return

    #...........................................................................
    # Methods
    @PYB11virtual
    def solve(self,
              func = "SolverFunction&",
              x = "std::vector<double>&"):
        """Solve the system of equations represented by function 'func' and initial guess 'x'.
Solution vector is returned in 'x'.
Returns the number of non-linear iterations taken."""
        return "size_t"

    #...........................................................................
    # Properties
    globalstrategy = PYB11property("int", "globalstrategy")
    fnormtol = PYB11property("double", "fnormtol")
    scsteptol = PYB11property("double", "scsteptol")
    numMaxIters = PYB11property("long int", "numMaxIters", "numMaxIters")
