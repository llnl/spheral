import inspect
from PYB11Generator import *
from FieldSpan import FieldSpan

#-------------------------------------------------------------------------------
# Add numeric operations to a FieldSpan
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("FieldSpan")
class ArithmeticFieldSpan:

    PYB11typedefs = """
  using SelfType = FieldSpan<%(Dimension)s, %(Value)s>;
  using Scalar = typename SelfType::Scalar;
  using ScalarFieldSpan = FieldSpan<%(Dimension)s, Scalar>;
"""

    def __iadd__(self):
        return

    def __isub__(self):
        return

    @PYB11pyname("__iadd__")
    def __iadd__V__(self, rhs="%(Value)s()"):
        return

    @PYB11pyname("__isub__")
    def __isub__V__(self, rhs="%(Value)s()"):
        return

    @PYB11implementation("[](SelfType& self, const ScalarFieldSpan& rhs) { return self *= rhs; }")
    @PYB11operator
    def __imul__(self, rhs="const ScalarFieldSpan&"):
        return

    @PYB11implementation("[](SelfType& self, const ScalarFieldSpan& rhs) { return self /= rhs; }")
    @PYB11operator
    def __itruediv__(self, rhs="const ScalarFieldSpan&"):
        return

    @PYB11pyname("__imul__")
    def __imul__S__(self, rhs="Scalar()"):
        return

    @PYB11pyname("__itruediv__")
    def __itruediv__S__(self, rhs="Scalar()"):
        return

    @PYB11const
    def localSumElements(self):
        "Return the sum of the elements in the FieldSpan local to each processor."
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
        "Enforce a floor on the values of the FieldSpan."
        return

    def applyMax(self):
        "Enforce a ceiling on the values of the FieldSpan."
        return

    @PYB11const
    def localMin(self):
        "Return the mimimum value in the FieldSpan local to each processor."
        return

    @PYB11const
    def localMax(self):
        "Return the maximum value in the FieldSpan local to each processor."
        return

#-------------------------------------------------------------------------------
# Inject base FieldSpan methods
#-------------------------------------------------------------------------------
PYB11inject(FieldSpan, ArithmeticFieldSpan)
