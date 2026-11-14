#-------------------------------------------------------------------------------
# TreeNeighbor Value object.  The inherited TreeNeighborView API is exposed
# directly here; TreeNeighborView itself is intentionally not a Python type.
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Neighbor import *

@PYB11template("Dimension")
class TreeNeighbor(Neighbor):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using NodeID = std::pair<size_t, size_t>;
    using LevelKey = typename TreeNeighbor<%(Dimension)s>::LevelKey;
    using CellKey = typename TreeNeighbor<%(Dimension)s>::CellKey;
    using Plane = GeomPlane<%(Dimension)s>;
"""

    def pyinit(self,
               dataBase="DataBase<%(Dimension)s>&",
               searchType=("const NeighborSearchType", "NeighborSearchType::GatherScatter"),
               kernelExtent=("const double", "2.0"),
               xmin=("const Vector&", "-10.0 * Vector::one()"),
               xmax=("const Vector&", " 10.0 * Vector::one()")):
        "Construct a TreeNeighbor attached to a DataBase."

    @PYB11virtual
    def resizeNodeExtents(self):
        return "void"

    @PYB11virtual
    def clear(self):
        return "void"

    @PYB11virtual
    def updateNodes(self):
        return "void"

    @PYB11virtual
    def reinitialize(self):
        return "void"

    @PYB11virtual
    @PYB11pycppname("reinitialize")
    def reinitialize1(self,
                      xmin="const Vector&",
                      xmax="const Vector&",
                      htarget="const Scalar"):
        return "void"

    @PYB11const
    def gridLevel(self, H="const SymTensor&"):
        return "LevelKey"

    @PYB11const
    @PYB11pycppname("gridLevel")
    def gridLevel1(self, h="const Scalar&"):
        return "LevelKey"

    @PYB11const
    def cellSize(self, levelID="const LevelKey"):
        return "Scalar"

    @PYB11const
    def nearestCellCenter(self, xi="const Vector&", Hi="const SymTensor&"):
        return "Vector"

    @PYB11const
    @PYB11pycppname("nearestCellCenter")
    def nearestCellCenter1(self, xi="const Vector&", hi="const Scalar&"):
        return "Vector"

    @PYB11const
    def occupied(self, xi="const Vector&", Hi="const SymTensor&"):
        return "bool"

    @PYB11const
    @PYB11pycppname("occupied")
    def occupied1(self, xi="const Vector&", hi="const Scalar&"):
        return "bool"

    @PYB11implementation("[](TreeNeighbor<%(Dimension)s>& self, const DataBase<%(Dimension)s>& dataBase) { self.setNodeExtents(dataBase); }")
    def setNodeExtents(self, dataBase="const DataBase<%(Dimension)s>&"):
        return "void"

    @PYB11implementation("[](TreeNeighbor<%(Dimension)s>& self, const DataBase<%(Dimension)s>& dataBase) { self.setInternalNodeExtents(dataBase); }")
    def setInternalNodeExtents(self, dataBase="const DataBase<%(Dimension)s>&"):
        return "void"

    @PYB11implementation("[](TreeNeighbor<%(Dimension)s>& self, const DataBase<%(Dimension)s>& dataBase) { self.setGhostNodeExtents(dataBase); }")
    def setGhostNodeExtents(self, dataBase="const DataBase<%(Dimension)s>&"):
        return "void"

    @PYB11implementation("""[](TreeNeighbor<%(Dimension)s>& self,
                                const DataBase<%(Dimension)s>& dataBase,
                                const std::vector<NodeID>& nodeIDs) {
                                  self.setNodeExtents(dataBase, nodeIDs);
                                }""")
    @PYB11pycppname("setNodeExtents")
    def setNodeExtents1(self,
                        dataBase="const DataBase<%(Dimension)s>&",
                        nodeIDs="const std::vector<NodeID>&"):
        return "void"

    @PYB11implementation("""[](const TreeNeighbor<%(Dimension)s>& self,
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

    @PYB11implementation("""[](const TreeNeighbor<%(Dimension)s>& self,
                                const Vector& position,
                                const SymTensor& H,
                                std::vector<NodeID>& masterList,
                                std::vector<NodeID>& coarseNeighbors,
                                const bool ghostConnectivity) {
                                  return self.setMasterList(position, H, masterList, coarseNeighbors, ghostConnectivity);
                                }""")
    def setMasterList(self,
                      position="const Vector&",
                      H="const SymTensor&",
                      masterList="std::vector<NodeID>&",
                      coarseNeighbors="std::vector<NodeID>&",
                      ghostConnectivity=("const bool", "false")):
        return "std::tuple<bool, size_t, size_t>"

    @PYB11implementation("""[](const TreeNeighbor<%(Dimension)s>& self,
                                const Vector& position,
                                std::vector<NodeID>& masterList,
                                std::vector<NodeID>& coarseNeighbors,
                                const bool ghostConnectivity) {
                                  return self.setMasterList(position, masterList, coarseNeighbors, ghostConnectivity);
                                }""")
    @PYB11pycppname("setMasterList")
    def setMasterList1(self,
                       position="const Vector&",
                       masterList="std::vector<NodeID>&",
                       coarseNeighbors="std::vector<NodeID>&",
                       ghostConnectivity=("const bool", "false")):
        return "std::tuple<bool, size_t, size_t>"

    @PYB11implementation("""[](const TreeNeighbor<%(Dimension)s>& self,
                                const Plane& enterPlane,
                                const Plane& exitPlane,
                                std::vector<NodeID>& masterList,
                                std::vector<NodeID>& coarseNeighbors) {
                                  return self.setMasterList(enterPlane, exitPlane, masterList, coarseNeighbors);
                                }""")
    @PYB11pycppname("setMasterList")
    def setMasterList2(self,
                       enterPlane="const Plane&",
                       exitPlane="const Plane&",
                       masterList="std::vector<NodeID>&",
                       coarseNeighbors="std::vector<NodeID>&"):
        return "std::tuple<bool, size_t, size_t>"

    @PYB11implementation("""[](const TreeNeighbor<%(Dimension)s>& self,
                                const Vector& position,
                                const SymTensor& H,
                                const size_t numCoarse,
                                const std::vector<NodeID>& coarseNeighbors,
                                std::vector<NodeID>& refineNeighbors) {
                                  if (refineNeighbors.size() < numCoarse) refineNeighbors.resize(numCoarse);
                                  return self.setRefineNeighborList(position, H, numCoarse, coarseNeighbors, refineNeighbors);
                                }""")
    def setRefineNeighborList(self,
                              position="const Vector&",
                              H="const SymTensor&",
                              numCoarse="const size_t",
                              coarseNeighbors="const std::vector<NodeID>&",
                              refineNeighbors="std::vector<NodeID>&"):
        return "size_t"

    @PYB11implementation("""[](const TreeNeighbor<%(Dimension)s>& self,
                                const Vector& position,
                                const size_t numCoarse,
                                const std::vector<NodeID>& coarseNeighbors,
                                std::vector<NodeID>& refineNeighbors) {
                                  if (refineNeighbors.size() < numCoarse) refineNeighbors.resize(numCoarse);
                                  return self.setRefineNeighborList(position, numCoarse, coarseNeighbors, refineNeighbors);
                                }""")
    @PYB11pycppname("setRefineNeighborList")
    def setRefineNeighborList1(self,
                               position="const Vector&",
                               numCoarse="const size_t",
                               coarseNeighbors="const std::vector<NodeID>&",
                               refineNeighbors="std::vector<NodeID>&"):
        return "size_t"

    @PYB11implementation("""[](const TreeNeighbor<%(Dimension)s>& self,
                                const LevelKey levelID,
                                const CellKey cellID,
                                std::vector<NodeID>& masterList,
                                std::vector<NodeID>& coarseNeighbors,
                                const bool ghostConnectivity) {
                                  return self.setTreeMasterList(levelID, cellID, masterList, coarseNeighbors, ghostConnectivity);
                                }""")
    def setTreeMasterList(self,
                          levelID="const LevelKey",
                          cellID="const CellKey",
                          masterList="std::vector<NodeID>&",
                          coarseNeighbors="std::vector<NodeID>&",
                          ghostConnectivity=("const bool", "false")):
        return "std::tuple<bool, size_t, size_t>"

    @PYB11static
    @PYB11implementation("[](const Scalar& H, const Scalar kernelExtent) { return TreeNeighbor<%(Dimension)s>::Hextent(H, kernelExtent); }")
    def Hextent(self, H="const Scalar&", kernelExtent="const Scalar"):
        return "Vector"

    @PYB11static
    @PYB11implementation("[](const SymTensor& H, const Scalar kernelExtent) { return TreeNeighbor<%(Dimension)s>::Hextent(H, kernelExtent); }")
    @PYB11pycppname("Hextent")
    def Hextent1(self, H="const SymTensor&", kernelExtent="const Scalar"):
        return "Vector"

    @PYB11const
    def dumpTree(self, globalTree="const bool"):
        return "std::string"

    @PYB11const
    def dumpTreeStatistics(self, globalTree="const bool"):
        return "std::string"

    @PYB11const
    def serialize(self, buffer="std::vector<char>&"):
        return "void"

    @PYB11implementation("""[](TreeNeighbor<%(Dimension)s>& self,
                                const std::vector<char>& buffer) {
                                  auto itr = buffer.begin();
                                  self.deserialize(itr, buffer.end());
                                }""")
    def deserialize(self, buffer="const std::vector<char>&"):
        return "void"

    @PYB11const
    def valid(self):
        return "bool"

    neighborSearchType = PYB11property("NeighborSearchType", "neighborSearchType",
                                       doc="The active neighbor search rule")
    kernelExtent = PYB11property("Scalar", "kernelExtent", "kernelExtent",
                                 doc="The interaction extent in eta space")
    xmin = PYB11property("const Vector&", "xmin", returnpolicy="reference_internal")
    xmax = PYB11property("const Vector&", "xmax", returnpolicy="reference_internal")
    boxLength = PYB11property("Scalar", "boxLength")
    numGridLevels = PYB11property("LevelKey", "numGridLevels")
    occupiedCells = PYB11property("std::vector<std::vector<CellKey>>", "occupiedCells")
