import inspect
from PYB11Generator import *
from FieldBase import *
from ArithmeticField import *

#-------------------------------------------------------------------------------
# Add min/max operations to a Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("Field")
class MinMaxField(FieldBase):

    PYB11typedefs = """
  using FieldType = Field<%(Dimension)s, %(Value)s>;
  using Scalar = typename FieldType::Scalar;
  using ScalarFieldType = Field<%(Dimension)s, Scalar>;
"""

    def applyScalarMin(self):
        "Enforce a double floor on the values of the Field."
        return

    def applyScalarMax(self):
        "Enforce a double ceiling on the values of the Field."
        return

#-------------------------------------------------------------------------------
# Inject base field methods
#-------------------------------------------------------------------------------
PYB11inject(ArithmeticField, MinMaxField)
