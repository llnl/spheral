import inspect
from PYB11Generator import *

#-------------------------------------------------------------------------------
# FieldListSet
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralFieldList")
class FieldListSet:

    # Constructors
    def pyinit(self,):
        "Default constructor"

    # Attributes
    ScalarFieldLists = PYB11readwrite(doc="The FieldList<Dim, double> set")
    VectorFieldLists = PYB11readwrite(doc="The FieldList<Dim, Vector> set")
    TensorFieldLists = PYB11readwrite(doc="The FieldList<Dim, Tensor> set")
    SymTensorFieldLists = PYB11readwrite(doc="The FieldList<Dim, SymTensor> set")
