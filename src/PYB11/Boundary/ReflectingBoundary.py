#-------------------------------------------------------------------------------
# ReflectingBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from PlanarBoundary import *
from BoundaryAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class ReflectingBoundary(PlanarBoundary):

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

    #............................................................................
    # Constructors
    def pyinit(self,
                dataBase = "const DataBase<%(Dimension)s>&",
                plane = "const Plane&"):
        "Construct with the reflecting plane"

    #............................................................................
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

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary10(self,
                             field = "Field<%(Dimension)s, RKCoefficients<%(Dimension)s>>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary21(self,
                             field = "Field<%(Dimension)s, std::vector<Vector>>&"):
        return "void"

    #............................................................................
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

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary10(self,
                          field = "Field<%(Dimension)s, RKCoefficients<%(Dimension)s>>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary21(self,
                          field = "Field<%(Dimension)s, std::vector<Vector>>&"):
        return "void"

    #...........................................................................
    @PYB11const
    def rkReflectOperator(self,
                          order = "const RKOrder",
                          useHessian = "const bool"):
        "The special reflection operator for RKCoefficients"
        return "const ReflectingBoundary<%(Dimension)s>::TransformationMatrix&"

    #...........................................................................
    # Properties
    reflectOperator = PYB11property(doc="The tensor reflection transformation")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, ReflectingBoundary)
PYB11inject(BoundaryAbstractMethods, ReflectingBoundary, virtual=True, pure_virtual=False)
