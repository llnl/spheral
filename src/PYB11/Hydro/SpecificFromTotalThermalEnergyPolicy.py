from PYB11Generator import *
from UpdatePolicyBase import *

@PYB11module("SpheralCompiledModules.SpheralHydro")
@PYB11template("Dimension")
class SpecificFromTotalThermalEnergyPolicy(UpdatePolicyBase):
    """SpecificFromTotalThermalEnergyPolicy -- An implementation of UpdatePolicyBase
specialized for the updating the specific thermal energy from the total
energy."""

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using KeyType = typename SpecificFromTotalThermalEnergyPolicy<%(Dimension)s>::KeyType;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self):
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

