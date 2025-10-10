#-------------------------------------------------------------------------------
# PeriodicBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from PlanarBoundary import *
from BoundaryAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class PeriodicBoundary(PlanarBoundary):

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
    def pyinit(self):
        "Default constructor"

    def pyinit1(self,
                plane1 = "const Plane&",
                plane2 = "const Plane&"):
        "Construct a periodic boundary mapping between the two (enter/exit) planes"

    #...........................................................................
    # Methods
    @PYB11virtual
    def cullGhostNodes(self,
                       flagSet = "const FieldList<%(Dimension)s, size_t>&",
                       old2newIndexMap = "FieldList<%(Dimension)s, size_t>&",
                       numNodesRemoved = "std::vector<size_t>&"):
        "Use a set of flags to cull out inactive ghost nodes."
        return "void"
    
    @PYB11virtual
    def reset(self,
              dataBase = "const DataBase<%(Dimension)s>&"):
        "Overridable hook for clearing out the boundary condition."
        return "void"

    @PYB11virtual
    @PYB11const
    def label(self):
        "Label for restart files"
        return "std::string"

    #............................................................................
    @PYB11pycppname("applyGhostBoundary")
    @PYB11const
    def applyGhostBoundary0(self,
                            fieldBase = "FieldBase<%(Dimension)s>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary9(self,
                            field = "Field<%(Dimension)s, FacetedVolume>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    #............................................................................
    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary9(self,
                         field = "Field<%(Dimension)s, FacetedVolume>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    #...........................................................................
    # Properties
    enter = PYB11property("const Plane&", "enterPlane", "setEnterPlane", doc="The first plane for periodic wrapping")
    exit = PYB11property("const Plane&", "exitPlane", "setExitPlane", doc="The second plane for periodic wrapping")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(BoundaryAbstractMethods, PeriodicBoundary, virtual=True, pure_virtual=False)
