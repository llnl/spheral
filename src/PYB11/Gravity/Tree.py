#-------------------------------------------------------------------------------
# Tree
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class Tree:

    "An implementation of a oct/quad-tree (3d/2d)"

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::Tensor Tensor;
  typedef typename %(Dimension)s::SymTensor SymTensor;
  typedef uint64_t CellKey;
  typedef uint32_t LevelKey;
  typedef std::pair<size_t, size_t> NodeID;
  typedef std::unordered_map<NodeID, std::vector<std::unordered_set<CellKey> > > CompletedCellSet;
"""

    #---------------------------------------------------------------------------
    # Cell nested type
    #---------------------------------------------------------------------------
    class Cell:

        PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::Tensor Tensor;
  typedef typename %(Dimension)s::SymTensor SymTensor;
  typedef uint64_t CellKey;
  typedef std::pair<size_t, size_t> NodeID;
"""

        # Constructors
        def pyinit(self):
            "default constructor"
            return

        def pyinit(self,
                   mi = "const double",
                   xi = "const Vector&",
                   vi = "const Vector&",
                   keyi = "const CellKey&"):
            "construct with mass, position, velocity, and key"
            return

        def pyinit(self,
                   mi = "const double",
                   xi = "const Vector&",
                   vi = "const Vector&",
                   keyi = "const CellKey&",
                   daughter = "const CellKey&"):
            "construct with mass, position, velocity, key, and key for daughter cell"
            return

        # Comparisons
        def __eq__(self):
            return
        def __lt__(self):
            return

        # attributes
        M = PYB11readwrite(doc="total mass")
        Mglobal = PYB11readwrite(doc="total mass across all processors")
        xcm = PYB11readwrite(doc="center of mass coordinate")
        vcm = PYB11readwrite(doc="center of mass velocity")
        rcm2cc2 = PYB11readwrite(doc="square of the distance between center of mass and geometric center")
        key = PYB11readwrite(doc="key for the cell")
        daughters = PYB11readwrite(doc="daughters of the cell", returnpolicy="reference")
        masses = PYB11readwrite(doc="masses of points that terminate in cell", returnpolicy="reference")
        positions = PYB11readwrite(doc="positions of points that terminate in cell", returnpolicy="reference")
        velocities = PYB11readwrite(doc="velocities of points that terminate in cell", returnpolicy="reference")
    #...........................................................................

    # Constructor
    def pyinit(self,
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "Construct with the given coordinate boundaries for the bounding box"
        return

    # Methods
    @PYB11const
    def cellCenter(self,
                   level = "const LevelKey&",
                   key = "const CellKey&"):
        "central position of cell"
        return "Vector"

    def lowerBound(self,
                   level = "const LevelKey&",
                   key = "const CellKey&"):
        "minimum coordinate in cell"
        return "Vector"

    def upperBound(self,
                   level = "const LevelKey&",
                   key = "const CellKey&"):
        "maximum coordinate in cell"
        return "Vector"


    def cellSize(self,
                 level = "const LevelKey&"):
        "extent of a cell on the given level"
        return "double"

    @PYB11implementation("""[](Tree<%(Dimension)s>& self, int i) -> py::list {
                              if (i < 0 or i > int(self.numLevels())) throw py::index_error();
                              py::list result;
                              const auto& cells = self[i];
                              for (const auto& x: cells) result.append(x.second);
                              return result;
                            }""")
    @PYB11const
    def __getitem__(self):
        "Return the cells on a grid level -- note this is a copy, so modifications are not reflected in the tree"
        return "py::list"

    @PYB11implementation("""[](Tree<%(Dimension)s>& self, const LevelKey level, const Vector& xi) -> py::tuple {
                                CellKey key, ix, iy, iz;
                                self.buildCellKey(level, xi, key, ix, iy, iz);
                                return py::make_tuple(key, ix, iy, iz);
                            }""")
    @PYB11const
    def buildCellKey(self,
                     level = "const LevelKey",
                     xi = "const Vector&"):
        "Return (cell_ key, x, iy, iz) for the given position on a grid level"
        return "py::tuple"

    # Properties
    num1dbits = PYB11readonly(static=True)
    max1dKey = PYB11readonly(static=True)
    xkeymask = PYB11readonly(static=True)
    ykeymask = PYB11readonly(static=True)
    zkeymask = PYB11readonly(static=True)
    xmin = PYB11property(doc="Minimum coordinate of bounding box for Tree")
    xmax = PYB11property(doc="Maximum coordinate of bounding box for Tree")
    numLevels = PYB11property("size_t", "numLevels", "numLevels", doc="number of levels in Tree")
