import inspect
from PYB11Generator import *

#-------------------------------------------------------------------------------
# FieldListSet
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralFieldList")
class FieldListSet:

    # Constructors
    def pyinit(self,):
        "Default constructor"

    # Attributes
    ScalarFieldLists = PYB11readwrite(doc="The FieldList<Dim, double> set", returnpolicy="reference_internal")
    VectorFieldLists = PYB11readwrite(doc="The FieldList<Dim, Vector> set", returnpolicy="reference_internal")
    TensorFieldLists = PYB11readwrite(doc="The FieldList<Dim, Tensor> set", returnpolicy="reference_internal")
    SymTensorFieldLists = PYB11readwrite(doc="The FieldList<Dim, SymTensor> set", returnpolicy="reference_internal")
