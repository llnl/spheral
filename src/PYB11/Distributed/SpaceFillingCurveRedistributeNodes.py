#-------------------------------------------------------------------------------
# SpaceFillingCurveRedistributeNodes
#-------------------------------------------------------------------------------
from PYB11Generator import *
from RedistributeNodes import *

@PYB11template("Dimension")
class SpaceFillingCurveRedistributeNodes(RedistributeNodes):
    """SpaceFillingCurveRedistributeNodes

This is an abstract base for the space filling curve family of 
repartitioners.  The assumption is that the descendent classes will provide
the computeHashedIndices method to assign unique keys to each point in the
order that that algorithm wants the points distributed."""

    PYB11typedefs = """
    typedef typename KeyTraits::Key Key;
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               minNodesPerDomainFraction = "const double",
               maxNodesPerDomainFraction = "const double",
               workBalance = ("const bool", "true"),
               localReorderOnly = ("const bool", "false")):
        "Constructor"

    #...........................................................................
    # Virtual methods
    @PYB11pure_virtual
    @PYB11const
    def computeHashedIndices(self,
                             dataBase = "const DataBase<%(Dimension)s>&"):
        "This is the required method for all descendant classes."
        return "FieldList<%(Dimension)s, Key> "

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
    def computeStepSize(self,
                        box = "const std::pair<Vector, Vector>&"):
        "Compute the cell size in each dimension."
        return "Vector"
  
    @PYB11const
    def buildIndex2IDPairs(self,
                           indices = "const FieldList<%(Dimension)s, Key>&",
                           domainNodes = "const std::vector<DomainNode<%(Dimension)s> >&"):
        """Stitch together the given indices and DomainNode list.
This returns the set sorted by the index."""
        return "std::vector<std::pair<Key, DomainNode<%(Dimension)s> > >"

    @PYB11const
    @PYB11implementation("""[](const SpaceFillingCurveRedistributeNodes<%(Dimension)s>& self,
                               const std::vector<Key>& indices,
                               const std::vector<int>& count,
                               const std::vector<Scalar>& work,
                               const Key lowerBound,
                               const Key maxUpperBound,
                               const Scalar workTarget,
                               const int minNodes,
                               const int maxNodes) {
                                   Key upperKey;
                                   int numNodes;
                                   auto key = self.findUpperKey(indices,
                                                                count,
                                                                work,
                                                                lowerBound,
                                                                maxUpperBound,
                                                                workTarget,
                                                                minNodes,
                                                                maxNodes,
                                                                upperKey,
                                                                numNodes);
                                   return py::make_tuple(key, upperKey, numNodes);
                               }""")
    def findUpperKey(self,
                     indices = "const std::vector<Key>&",
                     count = "const std::vector<int>&",
                     work = "const std::vector<Scalar>&",
                     lowerBound = "const Key",
                     maxUpperBound = "const Key",
                     workTarget = "const Scalar",
                     minNodes = "const int",
                     maxNodes = "const int"):
        "Find the hashed index the given amount of work above the specified lower bound."
        return "py::tuple"

    @PYB11const
    def numIndicesInRange(self,
                          indices = "const std::vector<Key>&",
                          count = "const std::vector<int>&",
                          lowerBound = "const Key",
                          upperBound = "const Key"):
        "Compute the global number of nodes in the given index range."
        return "int"

    @PYB11const
    def workInRange(self,
                    indices = "const std::vector<Key>&",
                    work = "const std::vector<Scalar>&",
                    lowerBound = "const Key",
                    upperBound = "const Key"):
        "Compute the global work for the nodes in the given index range."
        return "Scalar"

    @PYB11const
    @PYB11implementation("""[](const SpaceFillingCurveRedistributeNodes<%(Dimension)s>& self,
                               const std::vector<Key>& indices,
                               const std::vector<int>& count,
                               const std::vector<Scalar>& work,
                               const Key lowerBound,
                               const Key upperBound) {
                                   int countInRange;
                                   Scalar workInRange;
                                   self.workAndNodesInRange(indices, count, work, lowerBound, upperBound,
                                                            countInRange, workInRange);
                                   return py::make_tuple(countInRange, workInRange);
                               }""")
    def workAndNodesInRange(self,
                            indices = "const std::vector<Key>&",
                            count = "const std::vector<int>&",
                            work = "const std::vector<Scalar>&",
                            lowerBound = "const Key",
                            upperBound = "const Key"):
        "Combines numIndicesInRange and workInRange"
        return "py::tuple"

    @PYB11const
    def targetNumNodes(self,
                       numGlobal = "const int",
                       numProcs = "const int",
                       targetProc = "const int"):
        "Compute the number of nodes we want per process."
        return "int"

    @PYB11const
    def findNextIndex(self,
                      indices = "const std::vector<Key>&",
                      index = "const Key",
                      maxIndex = "const Key"):
        "Get the next (global) index following the given value."
        return "Key"

    @PYB11const
    def domainForIndex(self,
                       index = "const Key",
                       indexRanges = "const std::vector<std::pair<Key, Key> >&"):
        "Find the domain for the given index given the set of index ranges for each processor."
        return "int"

    #...........................................................................
    # Properties
    minNodesPerDomainFraction = PYB11property("double")
    maxNodesPerDomainFraction = PYB11property("double")
    workBalance = PYB11property("bool")
    localReorderOnly = PYB11property("bool")
