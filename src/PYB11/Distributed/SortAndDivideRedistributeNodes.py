#-------------------------------------------------------------------------------
# SortAndDivideRedistributeNodes
#-------------------------------------------------------------------------------
from PYB11Generator import *
from RedistributeNodes import *

@PYB11template("Dimension")
class SortAndDivideRedistributeNodes(RedistributeNodes):
    """SortAndDivideRedistributeNodes -- (Re)domain decompose the nodes by using 
a sort and divide algorithm.

This is a template base class -- the actual dimension dependent objects
are: 
  SortAndDivideRedistributeNodes1d
  SortAndDivideRedistributeNodes2d
  SortAndDivideRedistributeNodes3d"""

    PYB11typedefs = """
    typedef typename KeyTraits::Key Key;
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Constructor"

    #...........................................................................
    # Methods
    @PYB11const
    def sortByPositions(self,
                        domainNodes = "std::list<DomainNode<%(Dimension)s> >&",
                        positionIndex = "const int"):
        """Sort the given set of DomainNodes by a position component (where positionIndex
maps as
  0 ==> x,
  1 ==> y,
  2 ==> z),
and return the result as a list<DomainNode>."""
        return "void"

    @PYB11const
    def popFrontNodes(self,
                      sortedCandidateNodes = "std::list<DomainNode<%(Dimension)s> >&",
                      targetDomainWork = "const double",
                      positionIndex = "const int"):
        """Given a sorted list of domain nodes, pop nodes off of the front building
a return list until the requested work is met."""
        return "std::list<DomainNode<%(Dimension)s> >"

    @PYB11const
    def shapeTensor(self,
                    domainNodes = "const std::vector<DomainNode<%(Dimension)s> >&"):
        "Compute the appropriate shape tensor for a set of domain nodes."
        return "typename SymTensor::EigenStructType"

    @PYB11const
    def rotateIntoShapeTensorFrame(self,
                                   shapeTensor = "const typename SymTensor::EigenStructType&",
                                   domainNodes = "std::vector<DomainNode<%(Dimension)s> >&"):
        """Apply the necessary rotation to the positions of the domain nodes to transform into
the primary frame of the given shape tensor."""
        return "void"

    @PYB11const
    def reduceDomainNodes(self,
                          nodes = "const std::vector<DomainNode<%(Dimension)s> >&",
                          targetProc = "const int"):
        "Reduce a vector<DomainNode> to the given processor."
        return "std::vector<DomainNode<%(Dimension)s> >"

    @PYB11const
    def broadcastDomainNodes(self,
                             nodes = "const std::vector<DomainNode<%(Dimension)s> >&",
                             sendProc = "const int"):
        "Broadcast a vector<DomainNode> from the given processor."
        return "std::vector<DomainNode<%(Dimension)s> >"
