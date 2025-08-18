import inspect
from PYB11Generator import *
from FieldBase import *
from ArithmeticField import *
from MinMaxFieldSpan import MinMaxFieldSpan

#-------------------------------------------------------------------------------
# Add min/max operations to a Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("Field")
class MinMaxField(FieldBase,
                  MinMaxFieldSpan):

    PYB11typedefs = """
    using SelfType = Field<%(Dimension)s, %(Value)s>;
    using ViewType = typename SelfType::ViewType;
    using Scalar = typename SelfType::Scalar;
    using ScalarFieldType = Field<%(Dimension)s, Scalar>;
    using ScalarFieldSpan = FieldSpan<%(Dimension)s, Scalar>;
"""

    def applyScalarMin(self):
        "Enforce a double floor on the values of the Field."
        return

    def applyScalarMax(self):
        "Enforce a double ceiling on the values of the Field."
        return

    @PYB11const
    def min(self):
        "Return the mimimum value in the FieldSpan."
        return

    @PYB11const
    def max(self):
        "Return the maximum value in the FieldSpan."
        return

#-------------------------------------------------------------------------------
# Inject base field methods
#-------------------------------------------------------------------------------
PYB11inject(ArithmeticField, MinMaxField)
