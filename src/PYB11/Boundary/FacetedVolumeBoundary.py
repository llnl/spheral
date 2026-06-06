#-------------------------------------------------------------------------------
# ReflectingBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from PlanarBoundary import *
from BoundaryAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class FacetedVolumeBoundary(Boundary):

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
               poly = "const FacetedVolume&",
               interiorBoundary = "const bool",
               useGhosts = ("const bool", "false")):
        "Constructor"

    #............................................................................
    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary1(self,
                            field = "Field<%(Dimension)s, int>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary2(self,
                            field = "Field<%(Dimension)s, Scalar>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary3(self,
                            field = "Field<%(Dimension)s, Vector>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary4(self,
                            field = "Field<%(Dimension)s, Tensor>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary5(self,
                            field = "Field<%(Dimension)s, SymTensor>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary6(self,
                            field = "Field<%(Dimension)s, ThirdRankTensor>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary7(self,
                            field = "Field<%(Dimension)s, FourthRankTensor>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary8(self,
                            field = "Field<%(Dimension)s, FifthRankTensor>&"):
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
    def enforceBoundary1(self,
                         field = "Field<%(Dimension)s, int>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary2(self,
                         field = "Field<%(Dimension)s, Scalar>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary3(self,
                         field = "Field<%(Dimension)s, Vector>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary4(self,
                         field = "Field<%(Dimension)s, Tensor>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary5(self,
                         field = "Field<%(Dimension)s, SymTensor>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary6(self,
                         field = "Field<%(Dimension)s, ThirdRankTensor>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary7(self,
                         field = "Field<%(Dimension)s, FourthRankTensor>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary8(self,
                         field = "Field<%(Dimension)s, FifthRankTensor>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"
    
    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary9(self,
                         field = "Field<%(Dimension)s, FacetedVolume>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    #...........................................................................
    # Methods
    @PYB11virtual
    def reset(self,
              dataBase = "const DataBase<%(Dimension)s>&"):
        "Overridable hook for clearing out the boundary condition."
        return "void"

    @PYB11virtual
    def cullGhostNodes(self,
                       flagSet = "const FieldList<%(Dimension)s, size_t>&",
                       old2newIndexMap = "FieldList<%(Dimension)s, size_t>&",
                       numNodesRemoved = "std::vector<size_t>&"):
        "Use a set of flags to cull out inactive ghost nodes."
        return "void"

    @PYB11const
    def reflectOperator(self, 
                        facetID = "unsigned"):
        "Return the effective reflection operator for the given facet"
        return "const Tensor&"

    #...........................................................................
    # Properties
    polyVolume = PYB11property(doc="The faceted volume defining this boundary")
    interiorBoundary = PYB11property(doc="Flag if this Boundary is interior or exterior to the problem")
    useGhosts = PYB11property(doc="Flag whether we should build ghost points or not")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(BoundaryAbstractMethods, FacetedVolumeBoundary, virtual=True, pure_virtual=False)
