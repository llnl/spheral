from PYB11Generator import *

#-------------------------------------------------------------------------------
# FieldBase
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralField")
class FieldBase:

    # def pyinit(self, name="std::string"):
    #     "Construct with a name"

    # def pyinit1(self, name="std::string", nodeList="const NodeList<%(Dimension)s>&"):
    #     "Construct with a name and NodeList"

    @PYB11const
    def nodeList(self):
        "The NodeList this FieldBase is registered with."
        return "const NodeList<%(Dimension)s>&"

    @PYB11pure_virtual
    @PYB11const
    def size(self):
        "Number of elements"
        return "size_t"

    @PYB11pure_virtual
    def Zero(self):
        "Set all element values equal to zero"
        return "void"

    @PYB11pure_virtual
    def setNodeList(self, nodeList="const NodeList<%(Dimension)s>&"):
        "Register this Field with the given NodeList"
        return "void"

    @PYB11const
    def packValues(self, nodeIDs="const std::vector<size_t>&"):
        "Serialize the indicated elements into a vector<char>"
        return "std::vector<char>"

    @PYB11pure_virtual
    def unpackValues(self,
                     nodeIDs="const std::vector<size_t>&",
                     buffer = "const std::vector<char>&"):
        "Deserialize values from the given buffer"
        return "void"

    #...........................................................................
    # Properties
    name = PYB11property("std::string", getter="name", setter="name", doc="Name for the field")
