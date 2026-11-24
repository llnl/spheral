#-------------------------------------------------------------------------------
# VoronoiRedistributeNodes
#-------------------------------------------------------------------------------
from PYB11Generator import *
from RedistributeNodes import *

@PYB11template("Dimension")
class VoronoiRedistributeNodes(RedistributeNodes):
    """VoronoiRedistributeNodes

This algorithm uses the Voronoi tessellation to decide how to domain 
decompose our points.  The idea is to relax a set of generator points into
the SPH node distribution -- the generators are attracted to the SPH points
repelled by one and other.  These generator points then become the seeds to
draw the Voronoi tessellation about, each cell of which then represents a 
computational domain."""

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               workBalance = ("const bool", "false"),
               balanceGenerators = ("const bool", "true"),
               tolerance = ("const double", "1.0e-4"),
               maxIterations = ("const unsigned", "200")):
        "Constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def redistributeNodes(self,
                          dataBase = "DataBase<%(Dimension)s>&",
                          boundaries = ("std::vector<Boundary<%(Dimension)s>*>", "std::vector<Boundary<%(Dimension)s>*>()")):
        """Given a Spheral++ data base of NodeLists, repartition it among the processors.
This is the method required of all descendent classes."""
        return "void"

    #...........................................................................
    # Methods
    @PYB11const
    def computeCentroids(self,
                         nodes = "const std::vector<DomainNode<%(Dimension)s> >&",
                         generators = "std::vector<Vector>&"):
        "Given the set of DomainNodes with their domain assignments, compute the (work weighted) domain centroids."
        return "void"

    @PYB11const
    @PYB11implementation("""[](const VoronoiRedistributeNodes<%(Dimension)s>& self,
                               const std::vector<Vector>& generators,
                               const std::vector<int>& generatorFlags,
                               std::vector<double>& generatorWork,
                               std::vector<DomainNode<%(Dimension)s> >& nodes) {
                                   double minWork;
                                   double maxWork;
                                   unsigned minNodes;
                                   unsigned maxNodes;
                                   self.assignNodesToGenerators(generators, generatorFlags, generatorWork, nodes,
                                                                minWork, maxWork, minNodes, maxNodes);
                                   return py::make_tuple(minWork, maxWork, minNodes, maxNodes);
                               }""")
    def assignNodesToGenerators(self,
                                generators = "const std::vector<Vector>&",
                                generatorFlags = "const std::vector<int>&",
                                generatorWork = "std::vector<double>&",
                                nodes = "std::vector<DomainNode<%(Dimension)s> >&"):
        """Assign the nodes to the given generator positions, simultaneously computing the total generator work load.
Returns tuple (minWork, maxWork, minNodes, maxNodes)"""
        return "py::tuple"

    @PYB11const
    def cullGeneratorNodesByWork(self,
                                 generators = "const std::vector<Vector>&",
                                 generatorWork = "const std::vector<double>&",
                                 targetWork = "const double",
                                 generatorFlags = "std::vector<int>&",
                                 nodes = "std::vector<DomainNode<%(Dimension)s> >&"):
        "Look for any generator that has too much work, and unassign it's most distant nodes."
        return "void"

    @PYB11const
    def findNeighborGenerators(self,
                               igen = "const size_t",
                               generators = "const std::vector<Vector>&"):
        "Find the adjacent generators in the Voronoi diagram."
        return "std::vector<size_t>"

    #...........................................................................
    # Properties
    workBalance = PYB11property("bool")
    balanceGenerators = PYB11property("bool")
    tolerance = PYB11property("double")
    maxIterations = PYB11property("unsigned")
