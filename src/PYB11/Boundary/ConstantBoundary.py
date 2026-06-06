#-------------------------------------------------------------------------------
# ConstantBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from BoundaryAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class ConstantBoundary(Boundary):
    """ConstantBoundary -- Take a snapshot of the state of a set of nodes on 
NodeList, and create ghost nodes with that state from that time on.

This boundary is very specialized -- it explicitly works on only one 
NodeList.
"""

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename %(Dimension)s::FourthRankTensor FourthRankTensor;
    typedef typename %(Dimension)s::FifthRankTensor FifthRankTensor;
    typedef typename %(Dimension)s::FacetedVolume FacetedVolume;
    typedef GeomPlane<%(Dimension)s> Plane;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               dataBase = "DataBase<%(Dimension)s>&",
               nodeList = "NodeList<%(Dimension)s>&",
               nodeIndices = "const std::vector<size_t>&",
               denialPlane = "const Plane&"):
        "Construct a constant boundary for the specified nodes, including plane nodes are not allowed through"

    @PYB11implementation("""[](DataBase<%(Dimension)s>& db,
                               NodeList<%(Dimension)s>& nodes,
                               py::list nodeIndices,
                               const Plane& denialPlane) {
                                   return std::make_unique<ConstantBoundary<%(Dimension)s>>(db, nodes, Spheral::PYB11utils::from_list<size_t>(nodeIndices), denialPlane);
                               }""")
    def pyinit1(self,
                dataBase = "DataBase<%(Dimension)s>&",
                nodeList = "NodeList<%(Dimension)s>&",
                nodeIndices = "const std::vector<size_t>&",
                denialPlane = "const Plane&"):
        "Construct a constant boundary for the specified nodes, including plane nodes are not allowed through"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary(self,
                           fieldBase = "FieldBase<%(Dimension)s>&"):
        return "void"

    @PYB11virtual
    @PYB11const
    def enforceBoundary(self,
                        fieldBase = "FieldBase<%(Dimension)s>&"):
        return "void"

    @PYB11virtual
    def initializeProblemStartup(self,
                                 final = "const bool"):
        return "void"

    #...........................................................................
    # Properties
    nodeList = PYB11property("const NodeList<%(Dimension)s>&", "nodeList", doc="The NodeList this boundary applies to")
    nodeIndices = PYB11property("std::vector<size_t>", "nodeIndices", doc="The nodes this boundary is in control of")
    numConstantNodes = PYB11property("size_t", "numConstantNodes", doc="The number of nodes we are controlling")
    reflectOperator = PYB11property("const Tensor", "reflectOperator", doc="The tensor reflection transformation")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(BoundaryAbstractMethods, ConstantBoundary, virtual=True, pure_virtual=False)
PYB11inject(RestartMethods, ConstantBoundary)
