#-------------------------------------------------------------------------------
# ConnectivityMap
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
@PYB11holder("std::shared_ptr")
class ConnectivityMap:

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef NodeList<%(Dimension)s> NodeListType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    #...........................................................................
    # Methods
    def patchConnectivity(self,
                          flags = "const FieldList<%(Dimension)s, size_t>&",
                          old2new = "const FieldList<%(Dimension)s, size_t>&"):
        "Patch the connectivity information"
        return "void"

    def removeConnectivity(self,
                           neighborsToCut = "const FieldList<%(Dimension)s, std::vector<std::vector<int>>>&"):
        """Remove connectivity between neighbors.
Note this method assumes neighbor info is symmetric, and removes the pair connectivity for each
member of a pair (maintaining symmetry)."""
        return "void"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def connectivityForNode(self,
                            nodeList = "const NodeListType*",
                            nodeID = "const int"):
        "Get the set of neighbors for the given (internal!) node in the given NodeList."
        return "const std::vector<std::vector<int>>&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    @PYB11pycppname("connectivityForNode")
    def connectivityForNode1(self,
                             nodeListID = "const int",
                             nodeID = "const int"):
        "Get the set of neighbors for the given (internal!) node in the given NodeList."
        return "const std::vector<std::vector<int>>&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def intersectionConnectivity(self,
                                 pair = "const NodePairIdxType&"):
        """Get the pre-computed intersection connectivity for points if it was requested.
Note, this is different than what we expect for overlap connectivity: in this
method the intersection points (k) are all points that points (i,j) have in
common when (i,j) are ALSO neighbors.  Overlap connectivity may exist for
(i,j) even if (i,j) are not neighbors, and this set will miss such points."""
        return "const std::vector<std::vector<int>>&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def overlapConnectivityForNode(self,
                                   nodeList = "const NodeListType*",
                                   nodeID = "const int"):
        "The set of points that have non-zero overlap with the given point."
        return "const std::vector<std::vector<int>>&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    @PYB11pycppname("overlapConnectivityForNode")
    def overlapConnectivityForNode1(self,
                                    nodeListID = "const int",
                                    nodeID = "const int"):
        "The set of points that have non-zero overlap with the given point."
        return "const std::vector<std::vector<int>>&"

    @PYB11const
    def connectivityIntersectionForNodes(self,
                                         nodeListi = "const int",
                                         i = "const int",
                                         nodeListj = "const int",
                                         j = "const int",
                                         position = ("const FieldList<%(Dimension)s, Vector>&", "FieldList<%(Dimension)s, Vector>()")):
        "Compute the common neighbors for a pair of nodes."
        return "std::vector< std::vector<int> >"

    @PYB11const
    def connectivityUnionForNodes(self,
                                  nodeListi = "const int",
                                  i = "const int",
                                  nodeListj ="const int",
                                  j = "const int"):
        "Compute the union of neighbors for a pair of nodes."
        return "std::vector< std::vector<int> >"

    @PYB11const
    def numNeighborsForNode(self,
                            nodeListPtr = "const NodeListType*",
                            nodeID = "const int"):
        "Compute the number of neighbors for the given node."
        return "size_t"

    @PYB11const
    @PYB11pycppname("numNeighborsForNode")
    def numNeighborsForNode1(self,
                             nodeListID = "const int",
                             nodeID = "const int"):
        "Compute the number of neighbors for the given node."
        return "size_t"

    @PYB11const
    def globalConnectivity(self,
                           boundaries = "std::vector<Boundary<%(Dimension)s>*>&"):
        "Return the connectivity in terms of global node IDs."
        return "std::map<int, std::vector<int> >"

    @PYB11const
    def calculatePairInteraction(self,
                                 nodeListi = "const int",
                                 i = "const int", 
                                 nodeListj = "const int",
                                 j = "const int",
                                 firstGhostNodej ="const int"):
        "Function to determine if given node information (i and j), if the pair should already have been calculated by iterating over each others neighbors."
        return "bool"

    @PYB11const
    def numNodes(self, nodeList="const int"):
        "Return the number of nodes we should walk for the given NodeList."
        return "size_t"

    @PYB11const
    def ithNode(self,
                nodeList = "const int",
                index = "const int"):
        "The ith node (ordered) in the given NodeList."
        return "int"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def nodeList(self, index="const int"):
        "Get the ith NodeList or FluidNodeList."
        return "const NodeListType&"

    @PYB11const
    def nodeListIndex(self,  nodeList="const NodeListType*"):
        "Return which NodeList index in order the given one would be in our connectivity."
        return "unsigned"

    @PYB11const
    def valid(self):
        "Check that the internal data structure is valid."
        return "bool"

    #...........................................................................
    # Properties
    buildGhostConnectivity = PYB11property(doc="Are we building connectivity for ghost nodes?")
    buildOverlapConnectivity = PYB11property(doc="Are we building connectivity for nodes that overlap?")
    buildIntersectionConnectivity = PYB11property(doc="Are we building the connectivity intersection for nodes that interact?")
    nodeLists = PYB11property("const std::vector<const NodeListType*>&", "nodeLists",
                              returnpolicy="reference",
                              doc="The set of NodeLists we have connectivity for")
    nodePairList = PYB11property(returnpolicy="reference",
                                 doc="The connectivity as a set of (nodeListi, i, nodeListj, j)")
    coupling = PYB11property("const NodeCoupling&", returnpolicy="reference",
                             doc="The coupling functor for pairs of nodes")
