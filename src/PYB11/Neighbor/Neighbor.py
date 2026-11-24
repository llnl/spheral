#-------------------------------------------------------------------------------
# Neighbor -- CPU Value interface shared by Neighbor implementations.
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class Neighbor:

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using NodeID = std::pair<size_t, size_t>;
"""

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def dataBase(self):
        "The DataBase owning this Neighbor."
        return "const DataBase<%(Dimension)s>&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def nodeExtent(self):
        "The node extent FieldList maintained by this Neighbor."
        return "const FieldList<%(Dimension)s, Vector>&"

    @PYB11pure_virtual
    def resizeNodeExtents(self):
        return "void"

    @PYB11pure_virtual
    def clear(self):
        return "void"

    @PYB11pure_virtual
    def updateNodes(self):
        "Rebuild all Neighbor data."
        return "void"

    @PYB11implementation("""[](Neighbor<%(Dimension)s>& self,
                                std::vector<NodeID> nodeIDs) {
                                  self.updateNodes(std::span<NodeID>(nodeIDs));
                                }""")
    @PYB11pycppname("updateNodes")
    def updateNodes1(self, nodeIDs="std::vector<NodeID>"):
        "Rebuild Neighbor data for the supplied (NodeList,node) IDs."
        return "void"

    @PYB11pycppname("updateNodes")
    def updateNodes2(self, nodes="const NodeList<%(Dimension)s>&"):
        "Rebuild Neighbor data for one NodeList."
        return "void"

    @PYB11pure_virtual
    def reinitialize(self):
        "Reinitialize using the current bounds."
        return "void"

    @PYB11pure_virtual
    @PYB11pycppname("reinitialize")
    def reinitialize1(self,
                      xmin="const Vector&",
                      xmax="const Vector&",
                      htarget="const Scalar"):
        "Reinitialize using explicit bounds and target smoothing scale."
        return "void"
