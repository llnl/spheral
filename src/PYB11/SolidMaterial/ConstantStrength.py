#-------------------------------------------------------------------------------
# ConstantStrength
#-------------------------------------------------------------------------------
from PYB11Generator import *
from StrengthModel import *
from StrengthModelAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralSolidMaterial")
class ConstantStrength(StrengthModel):
    """ConstantStrength -- An implentation of StrengthModel returning constant
values for the shear modulus and yield strength."""

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ScalarField = Field<%(Dimension)s, Scalar>;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               mu0 = "const double",
               Y0 = "const double",
               muD = ("const double", "0.0"),
               YD = ("const double", "0.0")):
        "Construct with constant shear modulus (mu0) and yield strength (Y0), optionally specifying the damaged values as well (default to 0.0)."

    def pyinit1(self,
                mu0 = "const double",
                Y0 = "const double",
                eos = "const SolidEquationOfState<%(Dimension)s>&",
                muD = ("const double", "0.0"),
                YD = ("const double", "0.0")):
        """Construct with constant shear modulus (mu0), yield strength (Y0), and associated solid EOS.
Optionally you can specify damaged values (default to 0.0)."""

    #...........................................................................
    # Properties
    mu0 = PYB11property("double", doc="intact shear modulus")
    Y0 = PYB11property("double", doc="intact yield strength")
    muD = PYB11property("double", doc="damaged shear modulus")
    YD = PYB11property("double", doc="damaged yield strength")

#-------------------------------------------------------------------------------
# Inject Strength interface
#-------------------------------------------------------------------------------
PYB11inject(StrengthModelAbstractMethods, ConstantStrength, virtual=True, pure_virtual=False)
