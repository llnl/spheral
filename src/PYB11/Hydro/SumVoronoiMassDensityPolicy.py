from PYB11Generator import *
from UpdatePolicyBase import *

@PYB11module("SpheralCompiledModules.SpheralHydro")
@PYB11template("Dimension")
class SumVoronoiMassDensityPolicy(UpdatePolicyBase):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using KeyType = typename SumVoronoiMassDensityPolicy<%(Dimension)s>::KeyType;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self,
                W = "const TableKernel<%(Dimension)s>&",
                package = "const Physics<%(Dimension)s>&",
                rhoMin = "const double",
                rhoMax = "const double"):
        return

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def update(self,
               key = "const KeyType&",
               state = "State<%(Dimension)s>&",
               derivs = "StateDerivatives<%(Dimension)s>&",
               multiplier = "const double",
               t = "const double",
               dt = "const double"):
        "Update a FieldList assoicated with the given key"
        return "void"

    @PYB11virtual
    def updateAsIncrement(self,
                          key = "const KeyType&",
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&",
                          multiplier = "const double",
                          t = "const double",
                          dt = "const double"):
        """An alternate method to be called when you want to specify that the derivative information
should be assumed to not necessarily be properly time-centered, and therefore you should 
only use time advancement ideas, no "replace" or more sophisticated approaches.
Default to just calling the generic method."""
        return "void"
    
