import inspect
from PYB11Generator import *
from FieldView import FieldView

#-------------------------------------------------------------------------------
# Add numeric operations to a FieldView
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("FieldView")
class ArithmeticFieldView:

    PYB11typedefs = """
  using SelfType = FieldView<%(Dimension)s, %(Value)s>;
  using SelfFieldType = Field<%(Dimension)s, %(Value)s>;
  using Scalar = typename SelfType::Scalar;
  using ScalarFieldView = FieldView<%(Dimension)s, Scalar>;
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

    @PYB11implementation("[](SelfType& self, const ScalarFieldView& rhs) { return self *= rhs; }")
    @PYB11operator
    def __imul__(self, rhs="const ScalarFieldView&"):
        return

    @PYB11implementation("[](SelfType& self, const ScalarFieldView& rhs) { return self /= rhs; }")
    @PYB11operator
    def __itruediv__(self, rhs="const ScalarFieldView&"):
        return

    @PYB11pyname("__imul__")
    def __imul__S__(self, rhs="Scalar()"):
        return

    @PYB11pyname("__itruediv__")
    def __itruediv__S__(self, rhs="Scalar()"):
        return

    def applyMin(self):
        "Enforce a floor on the values of the FieldView."
        return

    def applyMax(self):
        "Enforce a ceiling on the values of the FieldView."
        return

    @PYB11const
    def localSumElements(self,
                          includeGhosts = ("bool", "false")):
        "Return the sum of the elements in the FieldView local to each processor."
        return "%(Value)s"

    @PYB11const
    def localMin(self,
                 includeGhosts = ("bool", "false")):
        "Return the mimimum value in the FieldView local to each processor."
        return "%(Value)s"

    @PYB11const
    def localMax(self,
                 includeGhosts = ("bool", "false")):
        "Return the maximum value in the FieldView local to each processor."
        return "%(Value)s"

#-------------------------------------------------------------------------------
# Inject base FieldView methods
#-------------------------------------------------------------------------------
PYB11inject(FieldView, ArithmeticFieldView)
