from PYB11Generator import *

#-------------------------------------------------------------------------------
# FieldBase
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralFieldList")
class FieldListBase:
    "Base class for FieldLists -- not much to implement in Python."

    # def pyinit(self, name="std::string"):
    #     "Construct with a name"

    # def pyinit1(self, name="std::string", nodeList="const NodeList<%(Dimension)s>&"):
    #     "Construct with a name and NodeList"
