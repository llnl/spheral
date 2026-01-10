#-------------------------------------------------------------------------------
# CollinsStrength
#-------------------------------------------------------------------------------
from PYB11Generator import *
from StrengthModel import *
from StrengthModelAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralSolidMaterial")
class CollinsStrength(StrengthModel):
    """CollinsStrength -- Implements a pressure dependent yield strength
model appropriate for geological materials.

Since this is solely a yield strength model it takes another StrengthModel
as an argument to compute the shear modulus.  Perhaps at some point we should
just split up the ideas of what provides shear modulus and yield strength?

   See Collins, Melosh, Ivanov, 2004 Appendix, MAPS"""

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ScalarField = Field<%(Dimension)s, Scalar>;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               shearModulusModel = "const StrengthModel<%(Dimension)s>&",
               mui = "const double",
               mud = "const double",
               Y0 = "const double",
               Ym = "const double"):
        """CollinsStrength constructor

mui : Coefficient of internal friction (intact material)
mud : Coefficient of internal friction (damaged material)
Y0  : Shear strength at zero pressure
yM  : von Mises plastic limit"""

    def pyinit1(self,
                shearModulusModel = "const StrengthModel<%(Dimension)s>&",
                mui = "const double",
                Y0 = "const double",
                Ym = "const double"):
        """CollinsStrength constructor

DEPRECATION WARNING: this version without mud is deprecated

mui : Coefficient of internal friction
Y0  : Shear strength at zero pressure
yM  : von Mises plastic limit"""

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def providesSoundSpeed(self):
        return "bool"

    @PYB11virtual
    @PYB11const
    def soundSpeed(self,
                   soundSpeed = "Field<%(Dimension)s, Scalar>&",
                   density = "const Field<%(Dimension)s, Scalar>&",
                   specificThermalEnergy = "const Field<%(Dimension)s, Scalar>&",
                   pressure = "const Field<%(Dimension)s, Scalar>&",
                   fluidSoundSpeed = "const Field<%(Dimension)s, Scalar>&",
                   damage = "const Field<%(Dimension)s, SymTensor>&"):
        return "void"

    #...........................................................................
    # Properties
    mui = PYB11property("double")
    mud = PYB11property("double")
    Y0 = PYB11property("double")
    Ym = PYB11property("double")

#-------------------------------------------------------------------------------
# Inject Strength interface
#-------------------------------------------------------------------------------
PYB11inject(StrengthModelAbstractMethods, CollinsStrength, virtual=True, pure_virtual=False)
