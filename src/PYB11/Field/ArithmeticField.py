import inspect
from PYB11Generator import *
from FieldBase import FieldBase
from Field import Field
from ArithmeticFieldView import ArithmeticFieldView

#-------------------------------------------------------------------------------
# Add numeric operations to a Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("Field")
class ArithmeticField(FieldBase,
                      ArithmeticFieldView):

    PYB11typedefs = """
    using SelfType = Field<%(Dimension)s, %(Value)s>;
    using ViewType = typename SelfType::ViewType;
    using Scalar = typename SelfType::Scalar;
    using ScalarFieldType = Field<%(Dimension)s, Scalar>;
    using ScalarFieldView = FieldView<%(Dimension)s, Scalar>;
"""

    def __add__(self):
        return

    def __sub__(self):
        return

    @PYB11pyname("__add__")
    def __add__V__(self, rhs="%(Value)s()"):
        return

    @PYB11pyname("__sub__")
    def __sub__V__(self, rhs="%(Value)s()"):
        return

    @PYB11implementation("[](const SelfType& self, const ScalarFieldType& rhs) { return self * rhs; }")
    @PYB11operator
    def __mul__(self, rhs="const ScalarFieldType&"):
        return "SelfType"

    @PYB11implementation("[](const SelfType& self, const ScalarFieldType& rhs) { return self / rhs; }")
    @PYB11operator
    def __truediv__(self, rhs="const ScalarFieldType&"):
        return "SelfType"

    @PYB11pyname("__mul__")
    def __mul__S__(self, rhs="Scalar()"):
        return

    @PYB11pyname("__truediv__")
    def __truediv__S__(self, rhs="Scalar()"):
        return

    @PYB11const
    def sumElements(self,
                    includeGhosts = ("bool", "false")):
        "Return the sum of the elements in the Field."
        return "%(Value)s"

    @PYB11const
    def min(self,
            includeGhosts = ("bool", "false")):
        "Return the mimimum value in the Field"
        return "%(Value)s"

    @PYB11const
    def max(self,
            includeGhosts = ("bool", "false")):
        "Return the maximum value in the Field"
        return "%(Value)s"

#-------------------------------------------------------------------------------
# Inject base field methods
#-------------------------------------------------------------------------------
PYB11inject(Field, ArithmeticField)
