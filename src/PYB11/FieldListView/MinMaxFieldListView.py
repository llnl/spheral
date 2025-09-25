from PYB11Generator import *
import FieldListView
from ArithmeticFieldListView import *

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

    def applyMin(self, rhs="const %(Value)s&"):
        "Enforce a %(Value)s floor on the values of the Field."
        return

    def applyMax(self, rhs="const %(Value)s&"):
        "Enforce a %(Value)s ceiling on the values of the Field."
        return

    def applyScalarMin(self, rhs="const Scalar"):
        "Enforce a double floor on the values of the Field."
        return

    def applyScalarMax(self, rhs="const Scalar"):
        "Enforce a double ceiling on the values of the Field."
        return

    @PYB11const
    def localMin(self,
                 includeGhosts = ("bool", "false")):
        "Return the mimimum value in the FieldListView local to each processor."
        return "%(Value)s"

    @PYB11const
    def localMax(self,
                 includeGhosts = ("bool", "false")):
        "Return the maximum value in the FieldListView local to each processor."
        return "%(Value)s"

#-------------------------------------------------------------------------------
# Inject FieldListView
#-------------------------------------------------------------------------------
PYB11inject(ArithmeticFieldListView, MinMaxFieldListView)
