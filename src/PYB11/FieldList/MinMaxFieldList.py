from PYB11Generator import *
from FieldList import FieldList
from FieldListBase import FieldListBase
from ArithmeticFieldList import ArithmeticFieldList

#-------------------------------------------------------------------------------
# Add min/max operations to a Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("FieldList")
class MinMaxFieldList(FieldListBase):

    PYB11typedefs = """
    using FieldListType = FieldList<%(Dimension)s, %(Value)s>;
    using FieldType = Field<%(Dimension)s, %(Value)s>;
    using NodeListType = NodeList<%(Dimension)s>;
    using Scalar = %(Dimension)s::Scalar;
    using Vector = %(Dimension)s::Vector;
    using SymTensor = %(Dimension)s::SymTensor;
    using ViewType = typename FieldListType::ViewType;
"""

    def applyMin(self, rhs="const %(Value)s&"):
        "Enforce a %(Value)s floor on the values of the Field."
        return

    def applyMax(self, rhs="const %(Value)s&"):
        "Enforce a %(Value)s ceiling on the values of the Field."
        return

    def applyScalarMin(self):
        "Enforce a double floor on the values of the Field."
        return

    def applyScalarMax(self):
        "Enforce a double ceiling on the values of the Field."
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

    @PYB11const
    def min(self,
            includeGhosts = ("bool", "false")):
        "Return the mimimum value in the Field."
        return "%(Value)s"

    @PYB11const
    def max(self,
            includeGhosts = ("bool", "false")):
        "Return the maximum value in the Field."
        return "%(Value)s"

#-------------------------------------------------------------------------------
# Inject FieldList
#-------------------------------------------------------------------------------
PYB11inject(ArithmeticFieldList, MinMaxFieldList)
