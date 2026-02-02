from PYB11Generator import *
from UpdatePolicyBase import *

@PYB11module("SpheralCompiledModules.SpheralDataBase")
@PYB11holder("std::shared_ptr")
@PYB11template("Dimension", "ValueType")
class FieldUpdatePolicy(UpdatePolicyBase):
    "FieldUpdatePolicy -- Base/interface class for the policies on how Field state variables are to be updated."

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def clonePerField(self):
        "Returns whether this policy should be cloned for each Field in a FieldList or not"
        return "bool"

    @PYB11virtual
    @PYB11const
    def serializeData(self,
                      buf = "std::vector<double>&",
                      key = "const KeyType&",
                      state = "const State<%(Dimension)s>&"):
        "Serialize the data in the Field to a buffer"
        return "void"

    @PYB11virtual
    @PYB11const
    def deserializeData(self,
                        buf = "const std::vector<double>&",
                        key = "const KeyType&",
                        state = "const State<%(Dimension)s>&",
                        offset = "const size_t"):
        "Deserialize the data in the Field from a buffer"
        return "size_t"
