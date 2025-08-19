import inspect
from PYB11Generator import *
from ArithmeticFieldView import *

#-------------------------------------------------------------------------------
# Add min/max operations to a FieldView
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("FieldView")
class MinMaxFieldView:

    PYB11typedefs = """
  using SelfType = FieldView<%(Dimension)s, %(Value)s>;
  using Scalar = typename SelfType::Scalar;
  using ScalarFieldView = FieldView<%(Dimension)s, Scalar>;
"""

    def applyScalarMin(self):
        "Enforce a double floor on the values of the FieldView."
        return

    def applyScalarMax(self):
        "Enforce a double ceiling on the values of the FieldView."
        return

#-------------------------------------------------------------------------------
# Inject base field methods
#-------------------------------------------------------------------------------
PYB11inject(ArithmeticFieldView, MinMaxFieldView)
