from PYB11Generator import *
from FieldList import FieldList as __FieldList
from FieldListBase import FieldListBase as __FieldListBase

#-------------------------------------------------------------------------------
# FieldList with numeric operations
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("FieldList")
class ArithmeticFieldList(__FieldListBase):

    PYB11typedefs = """
    using FieldListType = FieldList<%(Dimension)s, %(Value)s>;
    using FieldType = Field<%(Dimension)s, %(Value)s>;
    using NodeListType = NodeList<%(Dimension)s>;
    using Scalar = %(Dimension)s::Scalar;
    using Vector = %(Dimension)s::Vector;
    using SymTensor = %(Dimension)s::SymTensor;
    using ViewType = typename FieldListType::ViewType;
"""

    # @PYB11const
    # @PYB11cppname("operator()")
    # def valueat(self,
    #             position = "const %(Dimension)s::Vector&",
    #             W = "const TableKernel<%(Dimension)s>&"):
    #     "Return the interpolated value of the FieldList at a position."
    #     return "%(Value)s"

    def __add__(self):
        return

    def __sub__(self):
        return

    def __iadd__(self):
        return

    def __isub__(self):
        return

    @PYB11pyname("__add__")
    def __add__V(self, rhs="%(Value)s()"):
        return

    @PYB11pyname("__sub__")
    def __sub__V(self, rhs="%(Value)s()"):
        return

    @PYB11pyname("__iadd__")
    def __iadd__V(self, rhs="%(Value)s()"):
        return

    @PYB11pyname("__isub__")
    def __isub__V(self, rhs="%(Value)s()"):
        return

    def __imul__(self, rhs="double()"):
        return

    def __itruediv__(self, rhs="double()"):
        return

    @PYB11pyname("__imul__")
    def __imul__SFL(self, rhs="FieldList<%(Dimension)s, Scalar>()"):
        return

    @PYB11pyname("__itruediv__")
    def __itruediv__SFL(self, rhs="FieldList<%(Dimension)s, Scalar>()"):
        return

    @PYB11const
    def localSumElements(self,
                         includeGhosts = ("bool", "false")):
        "Return the sum of the elements in the FieldList local to each processor."
        return "%(Value)s"

    @PYB11const
    def sumElements(self,
                    includeGhosts = ("bool", "false")):
        "Return the sum of the elements in the Field."
        return "%(Value)s"

#-------------------------------------------------------------------------------
# Inject FieldList
#-------------------------------------------------------------------------------
PYB11inject(__FieldList, ArithmeticFieldList)
