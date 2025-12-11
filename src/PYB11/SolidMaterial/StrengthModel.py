#-------------------------------------------------------------------------------
# StrengthModel abstract class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from StrengthModelAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralSolidMaterial")
class StrengthModel:
    "Abstract base for strength models"

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ScalarField = Field<%(Dimension)s, Scalar>;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def providesSoundSpeed(self):
        return "bool"

    @PYB11virtual
    @PYB11const
    def providesBulkModulus(self):
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

    @PYB11virtual
    @PYB11const
    def bulkModulus(self,
                    bulkModulus = "Field<%(Dimension)s, Scalar>&",
                    massDensity = "const Field<%(Dimension)s, Scalar>&",
                    specificThermalEnergy = "const Field<%(Dimension)s, Scalar>&"):
        return "void"

    @PYB11virtual
    @PYB11const
    def meltSpecificEnergy(self,
                           meltSpecificEnergy = "Field<%(Dimension)s, Scalar>&",
                           density = "const Field<%(Dimension)s, Scalar>&",
                           specficThermalEnergy = "const Field<%(Dimension)s, Scalar>&"):
        return "void"

    @PYB11virtual
    @PYB11const
    def coldSpecificEnergy(self,
                           coldSpecificEnergy = "Field<%(Dimension)s, Scalar>&",
                           density = "const Field<%(Dimension)s, Scalar>&",
                           specficThermalEnergy = "const Field<%(Dimension)s, Scalar>&"):
        return "void"

#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(StrengthModelAbstractMethods, StrengthModel, virtual=False, pure_virtual=True)
