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
  using SelfFieldType = Field<%(Dimension)s, %(Value)s>;
  using Scalar = typename SelfType::Scalar;
  using ScalarFieldView = FieldView<%(Dimension)s, Scalar>;
"""

    #...........................................................................
    # Comparators
    def __gt__(self):
        return

    def __lt__(self):
        return

    def __ge__(self):
        return "bool"

    def __le__(self):
        return "bool"

    def __gt__(self, rhs="%(Value)s()"):
        "Greater than comparision with a %(Value)s"
        return "bool"

    def __lt__(self, rhs="%(Value)s()"):
        "Less than comparision with a %(Value)s"
        return "bool"

    def __ge__(self, rhs="%(Value)s()"):
        "Greater than or equal comparision with a %(Value)s"
        return "bool"

    def __le__(self, rhs="%(Value)s()"):
        "Less than or equal comparision with a %(Value)s"
        return "bool"

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
