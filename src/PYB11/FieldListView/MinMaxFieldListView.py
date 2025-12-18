from PYB11Generator import *
from ArithmeticFieldListView import ArithmeticFieldListView as __ArithmeticFieldListView

#-------------------------------------------------------------------------------
# Add min/max operations to a Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("FieldListView")
class MinMaxFieldListView:

    PYB11typedefs = """
    using FieldType = Field<%(Dimension)s, %(Value)s>;
    using FieldListType = FieldList<%(Dimension)s, %(Value)s>;
    using FieldViewType = FieldView<%(Dimension)s, %(Value)s>;
    using FieldListViewType = FieldListView<%(Dimension)s, %(Value)s>;
    using Scalar = %(Dimension)s::Scalar;
    using Vector = %(Dimension)s::Vector;
    using SymTensor = %(Dimension)s::SymTensor;
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

    def applyScalarMin(self, rhs="const Scalar"):
        "Enforce a double floor on the values of the Field."
        return

    def applyScalarMax(self, rhs="const Scalar"):
        "Enforce a double ceiling on the values of the Field."
        return

#-------------------------------------------------------------------------------
# Inject FieldListView
#-------------------------------------------------------------------------------
PYB11inject(__ArithmeticFieldListView, MinMaxFieldListView)
