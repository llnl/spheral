#-------------------------------------------------------------------------------
# NodeCoupling
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11module("SpheralCompiledModules.SpheralUtilities")
class NodeCoupling:
    "A functor base class encapsulating how we couple pairs of nodes."

    def pyinit(self):
        "Default constructor"

    @PYB11virtual
    @PYB11const
    @PYB11cppname("operator()")
    def __call__(self,
                 pair = "const NodePairIdxType&"):
        "Functional method to override for coupling (nodeListi, i) <-> (nodeListj, j)"
        return "double"

