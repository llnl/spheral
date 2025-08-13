import inspect
from PYB11Generator import *
from FieldBase import FieldBase
from Field import Field

#-------------------------------------------------------------------------------
# Add numeric operations to a Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("Field")
class ArithmeticField(FieldBase):

    PYB11typedefs = """
  using FieldType = Field<%(Dimension)s, %(Value)s>;
  using Scalar = typename FieldType::Scalar;
  using ScalarFieldType = Field<%(Dimension)s, Scalar>;
"""

    def __add__(self):
        return

    def __sub__(self):
        return

    def __iadd__(self):
        return

    def __isub__(self):
        return

    @PYB11pyname("__add__")
    def __add__V__(self, rhs="%(Value)s()"):
        return

    @PYB11pyname("__sub__")
    def __sub__V__(self, rhs="%(Value)s()"):
        return

    @PYB11pyname("__iadd__")
    def __iadd__V__(self, rhs="%(Value)s()"):
        return

    @PYB11pyname("__isub__")
    def __isub__V__(self, rhs="%(Value)s()"):
        return

    @PYB11implementation("[](const FieldType& self, const ScalarFieldType& rhs) { return self * rhs; }")
    @PYB11operator
    def __mul__(self, rhs="const ScalarFieldType&"):
        return "FieldType"

    @PYB11implementation("[](const FieldType& self, const ScalarFieldType& rhs) { return self / rhs; }")
    @PYB11operator
    def __truediv__(self, rhs="const ScalarFieldType&"):
        return "FieldType"

    @PYB11implementation("[](FieldType& self, const ScalarFieldType& rhs) { return self *= rhs; }")
    @PYB11operator
    def __imul__(self, rhs="const ScalarFieldType&"):
        return

    @PYB11implementation("[](FieldType& self, const ScalarFieldType& rhs) { return self /= rhs; }")
    @PYB11operator
    def __itruediv__(self, rhs="const ScalarFieldType&"):
        return

    @PYB11pyname("__mul__")
    def __mul__S__(self, rhs="Scalar()"):
        return

    @PYB11pyname("__truediv__")
    def __truediv__S__(self, rhs="Scalar()"):
        return

    @PYB11pyname("__imul__")
    def __imul__S__(self, rhs="Scalar()"):
        return

    @PYB11pyname("__itruediv__")
    def __itruediv__S__(self, rhs="Scalar()"):
        return

    @PYB11const
    def sumElements(self):
        "Return the sum of the elements in the Field."
        return

    @PYB11const
    def localSumElements(self):
        "Return the sum of the elements in the Field local to each processor."
        return

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

    def applyMin(self):
        "Enforce a floor on the values of the Field."
        return

    def applyMax(self):
        "Enforce a ceiling on the values of the Field."
        return

    @PYB11const
    def min(self):
        "Return the mimimum value in the Field."
        return

    @PYB11const
    def max(self):
        "Return the maximum value in the Field."
        return

    @PYB11const
    def localMin(self):
        "Return the mimimum value in the Field local to each processor."
        return

    @PYB11const
    def localMax(self):
        "Return the maximum value in the Field local to each processor."
        return

#-------------------------------------------------------------------------------
# Inject base field methods
#-------------------------------------------------------------------------------
PYB11inject(Field, ArithmeticField)
