from PYB11Generator import *
from UpdatePolicyBase import *

@PYB11module("SpheralHydro")
@PYB11template("Dimension")
class SpecificThermalEnergyPolicy(UpdatePolicyBase):
    """SpecificThermalEnergyPolicy -- An implementation of UpdatePolicyBase specialized
for the updating the specific thermal energy as a dependent quantity.

This version is specialized for the compatible energy discretization 
method described in
Owen, J. M. (2014). A compatibly differenced total energy conserving form of
SPH. International Journal for Numerical Methods in Fluids, 75(11), 749–774. """

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using KeyType = typename SpecificThermalEnergyPolicy<%(Dimension)s>::KeyType;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self,
                db = "const DataBase<%(Dimension)s>&"):
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
        """If the derivative stored values for the pair-accelerations has not been updated,
we need to just time advance normally."""
        return "void"

    @PYB11virtual
    @PYB11const
    def independent(self):
        return "bool"
