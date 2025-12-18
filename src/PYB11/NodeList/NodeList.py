from PYB11Generator import *
from RestartMethods import *

#-------------------------------------------------------------------------------
# NodeList template
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralNodeList")
@PYB11dynamic_attr
class NodeList:
    "Spheral NodeList base class in %(Dimension)s"

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef Field<%(Dimension)s, Scalar> ScalarField;
    typedef Field<%(Dimension)s, Vector> VectorField;
    typedef Field<%(Dimension)s, Tensor> TensorField;
    typedef Field<%(Dimension)s, SymTensor> SymTensorField;
"""

    def pyinit(self,
               name = "std::string",
               numInternal = ("size_t", "0u"),
               numGhost = ("size_t", "0u"),
               hmin = ("double", "1e-20"),
               hmax = ("double", "1e20"),
               hminratio = ("double", "0.1"),
               nPerh = ("double", "2.01"),
               maxNumNeighbors = ("size_t", "500u")):
        "Constructor for NodeList base class."
        return

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def mass(self):
        "The mass field"
        return "const ScalarField&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def positions(self):
        "The position field"
        return "const VectorField&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def velocity(self):
        "The velocity field"
        return "const VectorField&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def Hfield(self):
        "The H tensor field"
        return "const SymTensorField&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def work(self):
        "The CPU work field"
        return "ScalarField&"

    @PYB11pycppname("mass")
    def setmass(self, newValue="const ScalarField&"):
        "Set the mass field"
        return "void"

    @PYB11pycppname("positions")
    def setpositions(self, newValue="const VectorField&"):
        "Set the position field"
        return "void"

    @PYB11pycppname("velocity")
    def setvelocity(self, newValue="const VectorField&"):
        "Set the velocity field"
        return "void"

    @PYB11pycppname("Hfield")
    def setHfield(self, newValue="const SymTensorField&"):
        "Set the H tensor field"
        return "void"

    @PYB11pycppname("work")
    def setwork(self, newValue="const ScalarField&"):
        "Set the CPU work field"
        return "void"

    @PYB11const
    def Hinverse(self, result="SymTensorField&"):
        "Compute the inverse H field"
        return "void"

    @PYB11const
    def haveField(self, field="const FieldBase<%(Dimension)s>&"):
        "Test if the given field is defined on this NodeList"
        return "bool"

    @PYB11const
    def nodeType(self,
                 ID = "size_t"):
        "Return the classification of the given node"
        return "NodeType"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def neighbor(self):
        "Neighbor object associated with this NodeList"
        return "Neighbor<%(Dimension)s>&"

    def registerNeighbor(self, neighbor="Neighbor<%(Dimension)s>&"):
        "Associate a Neighbor object with this NodeList"
        return "void"

    def unregisterNeighbor(self):
        "Break the relation of this NodeList with it's Neighbor object"
        return "void"

    # Virtual methods
    @PYB11virtual
    def deleteNodes(self, nodeIDs="const std::vector<size_t>&"):
        "Delete the indicated nodes from this NodeList"
        return "void"

    @PYB11virtual
    def reorderNodes(self, newOrdering="const std::vector<size_t>&"):
        "Reorder the nodes to the given mapping"
        return "void"

    #...........................................................................
    # Comparison
    def __eq__(self):
        "Equivalence test with another NodeList"

    def __ne__(self):
        "Inequivalence test with another NodeList"

    @PYB11implementation("[](const NodeList<%(Dimension)s>& self) -> std::uintptr_t { return reinterpret_cast<std::uintptr_t>(&self); }")
    def __hash__(self):
        "Make NodeList objects hashable for Python"
        return "std::uintptr_t"

    #...........................................................................
    # Properties
    name = PYB11property("std::string", doc="Name of the NodeList")
    numNodes = PYB11property("size_t", doc="Total number of nodes in this NodeList")
    numInternalNodes = PYB11property("size_t", "numInternalNodes", "numInternalNodes", doc="Number of internal nodes in this NodeList")
    numGhostNodes = PYB11property("size_t", "numGhostNodes", "numGhostNodes", doc="Number of ghost nodes in this NodeList")
    numFields = PYB11property("size_t", doc="Number of fields defined on this NodeList")
    firstGhostNode = PYB11property("size_t", doc="Index of the first ghost node on this NodeList")
    nodesPerSmoothingScale = PYB11property("Scalar", "nodesPerSmoothingScale", "nodesPerSmoothingScale", 
                                           doc="The target number of nodes per smoothing scale")
    maxNumNeighbors = PYB11property("size_t", "maxNumNeighbors", "maxNumNeighbors",
                                    doc="The maximum number of neighbors per node allowed for this NodeList")
    hmin = PYB11property("Scalar", "hmin", "hmin", doc="Minimum allowed smoothing scale")
    hmax = PYB11property("Scalar", "hmax", "hmax", doc="Maximum allowed smoothing scale")
    hminratio = PYB11property("Scalar", "hminratio", "hminratio", 
                              doc="Minimum allowed ratio of min/max smoothing scale eigenvalues on each node")

#-------------------------------------------------------------------------------
# Inject the restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, NodeList)
