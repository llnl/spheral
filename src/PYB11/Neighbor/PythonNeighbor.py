#-------------------------------------------------------------------------------
# CPU-only Python-overridable Neighbor Value object.
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Neighbor import *

@PYB11template("Dimension")
class PythonNeighbor(Neighbor):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using NodeID = std::pair<size_t, size_t>;
    using Plane = GeomPlane<%(Dimension)s>;
"""

    def pyinit(self,
               dataBase="DataBase<%(Dimension)s>&",
               searchType=("const NeighborSearchType", "NeighborSearchType::GatherScatter"),
               kernelExtent=("const Scalar", "2.0")):
        "Construct a CPU-only Python-overridable Neighbor."

    @PYB11virtual
    def resizeNodeExtents(self):
        return "void"

    @PYB11virtual
    def clear(self):
        return "void"

    @PYB11pure_virtual
    def updateNodes(self):
        return "void"

    @PYB11pure_virtual
    def reinitialize(self):
        return "void"

    @PYB11pure_virtual
    @PYB11pycppname("reinitialize")
    def reinitialize1(self,
                      xmin="const Vector&",
                      xmax="const Vector&",
                      htarget="const Scalar"):
        return "void"

    @PYB11pure_virtual
    @PYB11const
    def setMasterList(self,
                      position="const Vector&",
                      H="const SymTensor&",
                      masterList="std::vector<NodeID>&",
                      coarseNeighbors="std::vector<NodeID>&",
                      ghostConnectivity=("const bool", "false")):
        return "std::tuple<bool, size_t, size_t>"

    @PYB11pure_virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList1(self,
                       position="const Vector&",
                       masterList="std::vector<NodeID>&",
                       coarseNeighbors="std::vector<NodeID>&",
                       ghostConnectivity=("const bool", "false")):
        return "std::tuple<bool, size_t, size_t>"

    @PYB11pure_virtual
    @PYB11const
    @PYB11pycppname("setMasterList")
    def setMasterList2(self,
                       enterPlane="const Plane&",
                       exitPlane="const Plane&",
                       masterList="std::vector<NodeID>&",
                       coarseNeighbors="std::vector<NodeID>&"):
        return "std::tuple<bool, size_t, size_t>"

    @PYB11pure_virtual
    @PYB11const
    def setRefineNeighborList(self,
                              position="const Vector&",
                              H="const SymTensor&",
                              numCoarse="const size_t",
                              coarseNeighbors="const std::vector<NodeID>&",
                              refineNeighbors="std::vector<NodeID>&"):
        return "size_t"

    @PYB11pure_virtual
    @PYB11const
    @PYB11pycppname("setRefineNeighborList")
    def setRefineNeighborList1(self,
                               position="const Vector&",
                               numCoarse="const size_t",
                               coarseNeighbors="const std::vector<NodeID>&",
                               refineNeighbors="std::vector<NodeID>&"):
        return "size_t"

    @PYB11implementation("[](PythonNeighbor<%(Dimension)s>& self, const DataBase<%(Dimension)s>& dataBase) { self.setNodeExtents(dataBase); }")
    def setNodeExtents(self, dataBase="const DataBase<%(Dimension)s>&"):
        return "void"

    @PYB11implementation("[](PythonNeighbor<%(Dimension)s>& self, const DataBase<%(Dimension)s>& dataBase) { self.setInternalNodeExtents(dataBase); }")
    def setInternalNodeExtents(self, dataBase="const DataBase<%(Dimension)s>&"):
        return "void"

    @PYB11implementation("[](PythonNeighbor<%(Dimension)s>& self, const DataBase<%(Dimension)s>& dataBase) { self.setGhostNodeExtents(dataBase); }")
    def setGhostNodeExtents(self, dataBase="const DataBase<%(Dimension)s>&"):
        return "void"

    @PYB11implementation("""[](PythonNeighbor<%(Dimension)s>& self,
                                const DataBase<%(Dimension)s>& dataBase,
                                const std::vector<NodeID>& nodeIDs) {
                                  self.setNodeExtents(dataBase, nodeIDs);
                                }""")
    @PYB11pycppname("setNodeExtents")
    def setNodeExtents1(self,
                        dataBase="const DataBase<%(Dimension)s>&",
                        nodeIDs="const std::vector<NodeID>&"):
        return "void"

    @PYB11implementation("""[](const PythonNeighbor<%(Dimension)s>& self,
                                FieldList<%(Dimension)s, Vector>& positions,
                                const Vector& minMasterPosition,
                                const Vector& maxMasterPosition,
                                const Vector& minMasterExtent,
                                const Vector& maxMasterExtent,
                                const size_t numCoarse,
                                std::vector<NodeID>& coarseList) {
                                  return self.precullList(positions.view(), minMasterPosition, maxMasterPosition,
                                                          minMasterExtent, maxMasterExtent, numCoarse, coarseList);
                                }""")
    def precullList(self,
                    positions="FieldList<%(Dimension)s, Vector>&",
                    minMasterPosition="const Vector&",
                    maxMasterPosition="const Vector&",
                    minMasterExtent="const Vector&",
                    maxMasterExtent="const Vector&",
                    numCoarse="const size_t",
                    coarseList="std::vector<NodeID>&"):
        return "size_t"

    @PYB11static
    @PYB11implementation("[](const Scalar& H, const Scalar kernelExtent) { return PythonNeighbor<%(Dimension)s>::Hextent(H, kernelExtent); }")
    def Hextent(self, H="const Scalar&", kernelExtent="const Scalar"):
        return "Vector"

    @PYB11static
    @PYB11implementation("[](const SymTensor& H, const Scalar kernelExtent) { return PythonNeighbor<%(Dimension)s>::Hextent(H, kernelExtent); }")
    @PYB11pycppname("Hextent")
    def Hextent1(self, H="const SymTensor&", kernelExtent="const Scalar"):
        return "Vector"

    neighborSearchType = PYB11property("NeighborSearchType", "neighborSearchType")
    kernelExtent = PYB11property("Scalar", "kernelExtent", "kernelExtent")
