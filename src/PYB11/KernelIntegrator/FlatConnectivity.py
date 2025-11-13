from PYB11Generator import *

@PYB11holder("std::shared_ptr")
@PYB11template("Dimension")
class FlatConnectivity:
    PYB11typedefs = """
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename std::array<int, %(Dimension)s::nDim> ArrayDim;
"""
    def pyinit(self):
        "Indexing for bilinear form"

    @PYB11const
    def indexingInitialized(self):
        return "bool"

    @PYB11const
    def overlapIndexingInitialized(self):
        return "bool"

    @PYB11const
    def globalIndexingInitialized(self):
        return "bool"
    
    @PYB11const
    def surfaceIndexingInitialized(self):
        return "bool"
    
    @PYB11const
    def boundaryInformationInitialized(self):
        return "bool"

    @PYB11const
    def firstGlobalIndex(self):
        return "int"

    @PYB11const
    def lastGlobalIndex(self):
        return "int"

    @PYB11const
    def numNodes(self):
        return "int"

    @PYB11const
    def numInternalNodes(self):
        return "int"

    @PYB11const
    def numGlobalNodes(self):
        return "int"

    @PYB11const
    def numBoundaryNodes(self):
        return "int"
    
    @PYB11const
    def nodeToLocal(self,
                    nodeListi = "const int",
                    nodei = "const int"):
        return "int"

    @PYB11const
    def localToNode(self,
                    locali = "const int"):
        return "std::pair<int, int>"

    @PYB11const
    def localToGlobal(self,
                      locali = "const int"):
        return "int"

    @PYB11const
    def numNeighbors(self,
                     locali = "const int"):
        return "int"

    @PYB11const
    def numOverlapNeighbors(self,
                            locali = "const int"):
        return "int"

    @PYB11const
    def numConstNeighbors(self,
                          locali = "const int"):
        return "int"

    @PYB11const
    def numConstOverlapNeighbors(self,
                                 locali = "const int"):
        return "int"

    @PYB11const
    def numNonConstNeighbors(self,
                             locali = "const int"):
        return "int"

    @PYB11const
    def numNonConstOverlapNeighbors(self,
                                    locali = "const int"):
        return "int"

    @PYB11const
    def localToFlat(self,
                    locali = "const int",
                    localj = "const int"):
        return "int"

    @PYB11const
    def localToFlatOverlap(self,
                           locali = "const int",
                           localj = "const int"):
        return "int"

    @PYB11const
    def flatToLocal(self,
                    locali = "const int",
                    flatj = "const int"):
        return "int"

    @PYB11const
    def flatOverlapToLocal(self,
                           locali = "const int",
                           flatj = "const int"):
        return "int"
    
    @PYB11const
    def isConstantBoundaryNode(self,
                               locali = "const int"):
        return "bool"

    @PYB11const
    @PYB11implementation("""[](const FlatConnectivity<%(Dimension)s>& self,
                               const int locali) {
                               std::vector<int> result;
                               self.neighborIndices(locali, result);
                               return PYB11utils::to_list(result); }""")
    def neighborIndices(self,
                        locali = "const int"):
        return "py::list"

    @PYB11const
    @PYB11implementation("""[](const FlatConnectivity<%(Dimension)s>& self,
                               const int locali) {
                               std::vector<int> result;
                               self.overlapNeighborIndices(locali, result);
                               return PYB11utils::to_list(result); }""")
    def overlapNeighborIndices(self,
                               locali = "const int"):
        return "py::list"

    @PYB11const
    @PYB11implementation("""[](const FlatConnectivity<%(Dimension)s>& self,
                               const int locali) {
                               std::vector<int> result;
                               self.constNeighborIndices(locali, result);
                               return PYB11utils::to_list(result); }""")
    def constNeighborIndices(self,
                             locali = "const int"):
        return "py::list"

    @PYB11const
    @PYB11implementation("""[](const FlatConnectivity<%(Dimension)s>& self,
                               const int locali) {
                               std::vector<int> result;
                               self.overlapConstNeighborIndices(locali, result);
                               return PYB11utils::to_list(result); }""")
    def overlapConstNeighborIndices(self,
                                    locali = "const int"):
        return "py::list"

    @PYB11const
    @PYB11implementation("""[](const FlatConnectivity<%(Dimension)s>& self,
                               const int locali) {
                               std::vector<int> result;
                               self.nonConstNeighborIndices(locali, result);
                               return PYB11utils::to_list(result); }""")
    def nonConstNeighborIndices(self,
                                locali = "const int"):
        return "PYB11utils::to_list(result)"

    @PYB11const
    @PYB11implementation("""[](const FlatConnectivity<%(Dimension)s>& self,
                               const int locali) {
                               std::vector<int> result;
                               self.overlapNonConstNeighborIndices(locali, result);
                               return PYB11utils::to_list(result); }""")
    def overlapNonConstNeighborIndices(self,
                                       locali = "const int"):
        return "py::list"

    @PYB11const
    @PYB11implementation("""[](const FlatConnectivity<%(Dimension)s>& self,
                               const int locali) {
                               std::vector<int> result;
                               self.globalNeighborIndices(locali, result);
                               return PYB11utils::to_list(result); }""")
    def globalNeighborIndices(self,
                              locali = "const int"):
        return "py::list"

    @PYB11const
    @PYB11implementation("""[](const FlatConnectivity<%(Dimension)s>& self,
                               const int locali) {
                               std::vector<int> result;
                               self.globalOverlapNeighborIndices(locali, result);
                               return PYB11utils::to_list(result); }""")
    def globalOverlapNeighborIndices(self,
                                     locali = "const int"):
        return "py::list"

    @PYB11const
    @PYB11implementation("""[](const FlatConnectivity<%(Dimension)s>& self,
                               const unsigned locali) {
                               std::vector<unsigned> v1, v2, v3, v4;
                               const auto res = self.uniqueNeighborIndices(locali, v1, v2, v3, v4);
                               return py::make_tuple(res, v1, v2, v3, v4); }""")
    def uniqueNeighborIndices(self,
                              locali = "const int"):
        return "py::tuple"
    
    @PYB11const
    def numSurfaces(self,
                    locali = "const int"):
        return "int"

    @PYB11const
    def surfaceIndex(self,
                     locali = "const int",
                     normal = "const Vector&"):
        return "int"

    @PYB11const
    def normal(self,
               locali = "const int",
               flats = "const int"):
        return "const Vector&"
    
    @PYB11const
    def numSurfacesForCell(self,
                           locali = "const int"):
        return "int"
    
    @PYB11const
    def surfaceIndexForCell(self,
                            locali = "const int",
                            flats = "const int"):
        return "int"

    def computeIndices(self,
                       dataBase = "const DataBase<%(Dimension)s>&"):
        return "void"

    def computeOverlapIndices(self,
                              dataBase = "const DataBase<%(Dimension)s>&"):
        return "void"
    
    def computeGlobalIndices(self,
                             dataBase = "const DataBase<%(Dimension)s>&",
                             boundaries = "const std::vector<Boundary<%(Dimension)s>*>&"):
        return "void"
    
    def computeSurfaceIndices(self,
                              dataBase = "const DataBase<%(Dimension)s>&",
                              state = "const State<%(Dimension)s>&"):
        return "void"

    def computeBoundaryInformation(self,
                                   dataBase = "const DataBase<%(Dimension)s>&",
                                   boundaries = "const std::vector<Boundary<%(Dimension)s>*>&"):
        return "void"
