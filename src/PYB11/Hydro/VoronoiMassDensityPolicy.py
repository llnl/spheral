from PYB11Generator import *
from UpdatePolicyBase import *

@PYB11module("SpheralCompiledModules.SpheralHydro")
@PYB11template("Dimension")
class VoronoiMassDensityPolicy(UpdatePolicyBase):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using KeyType = typename VoronoiMassDensityPolicy<%(Dimension)s>::KeyType;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self,
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

