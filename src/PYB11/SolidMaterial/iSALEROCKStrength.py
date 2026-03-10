#-------------------------------------------------------------------------------
# iSALEROCKStrength
#-------------------------------------------------------------------------------
from PYB11Generator import *
from StrengthModel import *
from StrengthModelAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralSolidMaterial")
class iSALEROCKStrength(StrengthModel):
    """iSALEROCKStrength -- Implements a pressure dependent yield strength
model appropriate for geological materials.

Since this is solely a yield strength model it takes another StrengthModel
as an argument to compute the shear modulus.  Perhaps at some point we should
just split up the ideas of what provides shear modulus and yield strength?

   See Collins, Melosh, Ivanov, 2004 Appendix, MAPS
       Raducan et al. 2020 (projectile shape effects paper)"""

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ScalarField = Field<%(Dimension)s, Scalar>;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               shearModulusModel = "const StrengthModel<%(Dimension)s>&",
               Yi0 = "const double",                                   # Intact strength at zero pressure     
               Yiinf = "const double",                                 # Intact strength at infinite pressure 
               fi = "const double",                                    # Intact internal friction coefficient 
               Yd0 = "const double",                                   # Damaged strength at zero pressure    
               Ydinf = "const double",                                 # Damaged strength at infinite pressure
               fd = "const double"):                                   # Damaged internal friction coefficient
        """iSALEROCKStrength constructor

Yi0    : Intact strength at zero pressure     
Yiinf  : Intact strength at infinite pressure 
fi,    : Intact internal friction coefficient 
Yd0    : Damaged strength at zero pressure    
Ydinf  : Damaged strength at infinite pressure
fd)    : Damaged internal friction coefficient"""

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
    Yi0 = PYB11property("double", doc="Intact strength at zero pressure")
    Yiinf = PYB11property("double", doc="Intact strength at infinite pressure")
    fi = PYB11property("double", doc="Intact internal friction coefficient")
    Yd0 = PYB11property("double", doc="Damaged strength at zero pressure")
    Ydinf = PYB11property("double", doc="Damaged strength at infinite pressure")
    fd = PYB11property("double", doc="Damaged internal friction coefficient")

#-------------------------------------------------------------------------------
# Inject Strength interface
#-------------------------------------------------------------------------------
PYB11inject(StrengthModelAbstractMethods, iSALEROCKStrength, virtual=True, pure_virtual=False)
