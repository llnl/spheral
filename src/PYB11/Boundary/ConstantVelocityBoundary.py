#-------------------------------------------------------------------------------
# ConstantVelocityBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from BoundaryAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class ConstantVelocityBoundary(Boundary):
    """ConstantVelocityBoundary -- A boundary condition to enforce a constant 
velocity on a given set of nodes.

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
               nodeList = "const NodeList<%(Dimension)s>&",
               nodeIndices = "const std::vector<size_t>&"):
        "Construct a constant velocity boundary for the specified nodes"

    @PYB11implementation("[](const NodeList<%(Dimension)s>& nodeList, py::list nodeIndices) { return std::make_unique<ConstantVelocityBoundary<%(Dimension)s>>(nodeList, Spheral::PYB11utils::from_list<size_t>(nodeIndices)); }")
    def pyinit1(self,
                nodeList = "const NodeList<%(Dimension)s>&",
                nodeIndices = "py::list"):
        "Construct a constant velocity boundary for the specified nodes"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def valid(self):
        return "bool"

    #...........................................................................
    # Properties
    nodeList = PYB11property("const NodeList<%(Dimension)s>&", "nodeList", doc="The NodeList this boundary applies to")
    nodeIndices = PYB11property("std::vector<size_t>", "nodeIndices", doc="The nodes this boundary is in control of")
    velocityCondition = PYB11property("std::vector<Vector>", "velocityCondition", doc="The velocities for the nodes we control")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(BoundaryAbstractMethods, ConstantVelocityBoundary, virtual=True, pure_virtual=False)
PYB11inject(RestartMethods, ConstantVelocityBoundary)
