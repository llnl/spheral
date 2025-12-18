#-------------------------------------------------------------------------------
# DistributedBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *

from Boundary import *
from BoundaryAbstractMethods import *

@PYB11template("Dimension")
class DistributedBoundary(Boundary):
    """DistributedBoundary -- Base class for distributed parallel boundary
conditions, connecting NodeLists across parallel domains."""


    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename %(Dimension)s::FourthRankTensor FourthRankTensor;
    typedef typename %(Dimension)s::FifthRankTensor FifthRankTensor;
    typedef typename %(Dimension)s::FacetedVolume FacetedVolume;

    typedef typename DistributedBoundary<%(Dimension)s>::DomainBoundaryNodes DomainBoundaryNodes;
    typedef std::map<int, DomainBoundaryNodes> DomainBoundaryNodeMap;
    typedef std::map<NodeList<%(Dimension)s>*, DomainBoundaryNodeMap> NodeListDomainBoundaryNodeMap;
"""

    #...........................................................................
    class DomainBoundaryNodes:
        sendNodes = PYB11readwrite()
        receiveNodes = PYB11readwrite()

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    #...........................................................................
    # Methods
    @PYB11const
    def communicatedNodeList(self,
                             nodeList = "const NodeList<%(Dimension)s>&"):
        "Test if the given NodeList is communicated on this domain or not."
        return "bool"

    @PYB11const
    def nodeListSharedWithDomain(self,
                                 nodeList = "const NodeList<%(Dimension)s>&",
                                 neighborDomainID = "int"):
        "Test if the given NodeList is communicated with the given domain."
        return "bool"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def domainBoundaryNodeMap(self,
                              nodeList = "const NodeList<%(Dimension)s>&"):
        return "const DomainBoundaryNodeMap&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def domainBoundaryNodes(self,
                            nodeList = "const NodeList<%(Dimension)s>&",
                            neighborDomainID = "int"):
        return "const DomainBoundaryNodes&"

    @PYB11const
    def communicatedProcs(self,
                          sendProcs = "std::vector<int>&",
			  recvProcs = "std::vector<int>&"):
        "Extract the current set of processors we're communicating with."
        return "void"

    #...........................................................................
    # Virtual methods
    @PYB11pure_virtual
    def setAllGhostNodes(self,
                         dataBase = "DataBase<%(Dimension)s>&"):
        "Descendent Distributed Neighbors are required to provide the setGhostNodes method for DataBases."
        return "void"

    @PYB11virtual
    def cullGhostNodes(self,
                       flagSet = "const FieldList<%(Dimension)s, size_t>&",
                       old2newIndexMap = "FieldList<%(Dimension)s, size_t>&",
                       numNodesRemoved = "std::vector<size_t>&"):
        "Override the Boundary method for culling ghost nodes."
        return "void"

    @PYB11virtual
    @PYB11const
    def finalizeGhostBoundary(self):
        "Override the base method to finalize ghost boundaries."
        return "void"

    @PYB11virtual
    @PYB11const
    def meshGhostNodes(self):
        "We do not want to use the parallel ghost nodes as generators."
        return "bool"

    #...........................................................................
    # Non-blocking exchanges
    @PYB11const
    def beginExchangeFieldFixedSize(self,
                                    field = "FieldBase<%(Dimension)s>&"):
        "Start a non-blocking Field exchange"
        return "void"

    @PYB11const
    def beginExchangeFieldVariableSize(self,
                                       field = "FieldBase<%(Dimension)s>&"):
        "Start a non-blocking Field exchange"
        return "void"

    def finalizeExchanges(self):
        "Force the exchanges which have been registered to execute."
        return "void"

    @PYB11const
    def unpackField(self,
                    field = "FieldBase<%(Dimension)s>&",
                    packedValues = "const std::list< std::vector<char> >&"):
        "Unpack a packed set of Field values back into the Field."
        return "void"

    def setControlAndGhostNodes(self):
        "Update the control and ghost nodes of the base class"
        return "void"

    #...........................................................................
    # Protected methods
    @PYB11protected
    @PYB11returnpolicy("reference_internal")
    def accessNodeListDomainBoundaryNodeMap(self):
        "Descendent classes get read/write access to the communication maps."
        return "NodeListDomainBoundaryNodeMap&"

    @PYB11protected
    @PYB11returnpolicy("reference_internal")
    def accessDomainBoundaryNodeMap(self,
                                    nodeList = "const NodeList<%(Dimension)s>&"):
        "Descendent classes get read/write access to the communication maps."
        return "DomainBoundaryNodeMap&"

    @PYB11protected
    @PYB11returnpolicy("reference_internal")
    def accessDomainBoundaryNodes(self,
                                  nodeList = "const NodeList<%(Dimension)s>&",
                                  neighborDomainID = "int"):
        return "DomainBoundaryNodes&"

    @PYB11protected
    @PYB11returnpolicy("reference_internal")
    def openDomainBoundaryNodes(self,
                                nodeListPtr = "NodeList<%(Dimension)s>*",
                                domainID = "const int"):
        """Convenience method to return an iterator to the DomainBoundaryNodes for the
given NodeList and domain pair.  If there isn't an entry for this pair already,
this method creates one."""
        return "DomainBoundaryNodes&"

    @PYB11protected
    def removeDomainBoundaryNodes(self,
                                  nodeListPtr = "NodeList<%(Dimension)s>*",
                                  domainID = "const int"):
        "Inverse of above -- remove the DomainBoundaryNodes for a NodeList/procID pair."
        return "void"

    @PYB11virtual
    @PYB11protected
    def reset(self,
              dataBase = "const DataBase<%(Dimension)s>&"):
        "Override the Boundary method for clearing the maps."
        return "void"

    @PYB11protected
    def buildReceiveAndGhostNodes(self,
                                  dataBase = "const DataBase<%(Dimension)s>&"):
        """This handy helper method will build receive and ghost nodes on all each
domain based on send nodes that have already been filled in."""
        return "void"

    #...........................................................................
    # Properties
    domainID = PYB11property("int")
    numDomains = PYB11property("int")
    nodeListDomainBoundaryNodeMap = PYB11property("const NodeListDomainBoundaryNodeMap&", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(BoundaryAbstractMethods, DistributedBoundary, virtual=True, pure_virtual=False)
