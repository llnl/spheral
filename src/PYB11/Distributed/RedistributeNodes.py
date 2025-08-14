#-------------------------------------------------------------------------------
# RedistributeNodes
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class RedistributeNodes:
    """RedistributeNodes -- An abstract base class for methods that repartition
the Spheral++ NodeLists among domains."""

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    #...........................................................................
    # Virtual methods
    @PYB11pure_virtual
    def redistributeNodes(self,
                          dataBase = "DataBase<%(Dimension)s>&",
                          boundaries = ("std::vector<Boundary<%(Dimension)s>*>", "std::vector<Boundary<%(Dimension)s>*>()")):
        """Given a Spheral++ data base of NodeLists, repartition it among the processors.
This is the method required of all descendent classes."""
        return "void"

    #...........................................................................
    # Methods
    @PYB11const
    def numGlobalNodes(self,
                       dataBase = "const DataBase<%(Dimension)s>&"):
        "Global number of nodes in DataBase."
        return "int"

    @PYB11pycppname("currentDomainDecomposition")
    @PYB11const
    def currentDomainDecomposition1(self,
                                    dataBase = "const DataBase<%(Dimension)s>&",
                                    globalNodeIDs = "const FieldList<%(Dimension)s, size_t>&"):
        "Calculate the current domain decomposition, and return it as a set of DomainNode identifiers."
        return "std::vector<DomainNode<%(Dimension)s> >"

    @PYB11pycppname("currentDomainDecomposition")
    @PYB11const
    def currentDomainDecomposition2(self,
                                    dataBase = "const DataBase<%(Dimension)s>&",
                                    globalNodeIDs = "const FieldList<%(Dimension)s, size_t>&",
                                    workPerNode = "const FieldList<%(Dimension)s, Scalar>&"):
        "Same as currentDomainDecomposition, but fills in work field in the DomainNodes."
        return "std::vector<DomainNode<%(Dimension)s> >"

    @PYB11const
    def enforceDomainDecomposition(self,
                                   nodeDistribution = "const std::vector<DomainNode<%(Dimension)s> >&",
                                   dataBase = "DataBase<%(Dimension)s>&"):
        "Given a desired domain decomposition (as a vector<DomainNode>), reassign nodes appropriately."
        return "void"

    @PYB11const
    def validDomainDecomposition(self,
                                 nodeDistribution = "const std::vector<DomainNode<%(Dimension)s> >&",
                                 dataBase = "const DataBase<%(Dimension)s>&"):
        """Test that the given domain decomposition is valid (all nodes accounted 
for, once and only once, etc.)."""
        return "bool"

    @PYB11const
    def workPerNode(self,
                    dataBase = "const DataBase<%(Dimension)s>&",
                    Hextent = "const double"):
        "Compute the work per node."
        return "FieldList<%(Dimension)s, Scalar>"

    @PYB11const
    def gatherDomainDistributionStatistics(self,
                                           work = "const FieldList<%(Dimension)s, Scalar>&"):
        """Gather up and print the statistics of the current domain distribution based on
a computed work field."""
        return "std::string"

    #...........................................................................
    # Protected methods
    @PYB11protected
    @PYB11const
    def packDomainNodes(self,
                        distribution = "const std::vector<DomainNode<%(Dimension)s> >&"):
        "Pack/unpack a vector<DomainNode> as a vector<double>, for use with MPI."
        return "std::vector<char>"

    @PYB11protected
    @PYB11const
    def unpackDomainNodes(self,
                          buf = "const std::vector<char>&"):
        "Pack/unpack a vector<DomainNode> as a vector<double>, for use with MPI."
        return "std::vector<DomainNode<%(Dimension)s> >"

    #...........................................................................
    # Properties
    domainID = PYB11property("int")
    numDomains = PYB11property("int")
    computeWork = PYB11property("bool")
