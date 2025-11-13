from PYB11Generator import *
from FieldListView import FieldListView as __FieldListView  # Prevent importing into main module namespace

#-------------------------------------------------------------------------------
# FieldListView with numeric operations
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("FieldListView")
@PYB11module("SpheralFieldListView")
class ArithmeticFieldListView:

    PYB11typedefs = """
    using FieldType = Field<%(Dimension)s, %(Value)s>;
    using FieldListType = FieldList<%(Dimension)s, %(Value)s>;
    using FieldViewType = FieldView<%(Dimension)s, %(Value)s>;
    using FieldListViewType = FieldListView<%(Dimension)s, %(Value)s>;
    using Scalar = %(Dimension)s::Scalar;
    using Vector = %(Dimension)s::Vector;
    using SymTensor = %(Dimension)s::SymTensor;
"""

    def __iadd__(self):
        return

    def __isub__(self):
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

    # @PYB11pyname("__imul__")
    # def __imul__SFL(self, rhs="const FieldListView<%(Dimension)s, Scalar>&"):
    #     return

    # @PYB11pyname("__itruediv__")
    # def __itruediv__SFL(self, rhs="const FieldListView<%(Dimension)s, Scalar>&"):
    #     return

    @PYB11const
    def localSumElements(self,
                         includeGhosts = ("bool", "false")):
        "Return the sum of the elements in the FieldListView local to each processor."
        return "%(Value)s"

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

    def applyMin(self, rhs="const %(Value)s&"):
        "Enforce a %(Value)s floor on the values of the FieldList."
        return

    def applyMax(self, rhs="const %(Value)s&"):
        "Enforce a %(Value)s ceiling on the values of the FieldList."
        return

    @PYB11const
    def localMin(self,
                 includeGhosts = ("bool", "false")):
        "Return the mimimum value in the FieldList local to each processor."
        return "%(Value)s"

    @PYB11const
    def localMax(self,
                 includeGhosts = ("bool", "false")):
        "Return the maximum value in the FieldList local to each processor."
        return "%(Value)s"

#-------------------------------------------------------------------------------
# Inject FieldListView
#-------------------------------------------------------------------------------
PYB11inject(__FieldListView, ArithmeticFieldListView)
