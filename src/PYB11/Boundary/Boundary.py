#-------------------------------------------------------------------------------
# Boundary base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from BoundaryAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralBoundary")
class Boundary:

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename %(Dimension)s::FourthRankTensor FourthRankTensor;
    typedef typename %(Dimension)s::FifthRankTensor FifthRankTensor;
    typedef typename %(Dimension)s::FacetedVolume FacetedVolume;
    typedef typename Boundary<%(Dimension)s>::BoundaryNodes BoundaryNodes;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    #...........................................................................
    # Virtual methods
    @PYB11pycppname("applyGhostBoundary")
    @PYB11const
    def applyGhostBoundary0(self,
                            fieldBase = "FieldBase<%(Dimension)s>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundaryInt(self,
                              field = "Field<%(Dimension)s, int>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundarySizeT(self,
                                field = "Field<%(Dimension)s, size_t>&"):
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

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary10(self,
                            field = "Field<%(Dimension)s, std::vector<Scalar>>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"
    
    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary11(self,
                            field = "Field<%(Dimension)s, std::vector<Vector>>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"
    
    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary12(self,
                             field = "Field<%(Dimension)s, RKCoefficients<%(Dimension)s>>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    #............................................................................
    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary0(self,
                         fieldBase = "FieldBase<%(Dimension)s>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundaryInt(self,
                           field = "Field<%(Dimension)s, int>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundarySizeT(self,
                             field = "Field<%(Dimension)s, size_t>&"):
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

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary10(self,
                         field = "Field<%(Dimension)s, std::vector<Scalar>>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"
    
    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary11(self,
                         field = "Field<%(Dimension)s, std::vector<Vector>>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"
    
    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary12(self,
                          field = "Field<%(Dimension)s, RKCoefficients<%(Dimension)s>>&"):
        "Apply the boundary condition to the violation node values in the given Field."
        return "void"

    #............................................................................
    @PYB11virtual
    def setAllGhostNodes(self,
                         dataBase = "DataBase<%(Dimension)s>&"):
        "Recreate ghost nodes for this Boundary for all NodeLists in the DataBase."
        return "void"

    @PYB11virtual
    def setAllViolationNodes(self,
                             dataBase = "DataBase<%(Dimension)s>&"):
        "Select any nodes based in the FluidNodeLists in the given DataBase that are in violation of boundary condition."
        return "void"

    @PYB11virtual
    def cullGhostNodes(self,
                       flagSet = "const FieldList<%(Dimension)s, size_t>&",
                       old2newIndexMap = "FieldList<%(Dimension)s, size_t>&",
                       numNodesRemoved = "std::vector<size_t>&"):
        "Use a set of flags to cull out inactive ghost nodes."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary20(self,
                             field = "Field<%(Dimension)s, std::vector<Scalar>>&"):
        "We provide default copies for arrays of values, but descendants can override these."
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary21(self,
                             field = "Field<%(Dimension)s, std::vector<Vector>>&"):
        "We provide default copies for arrays of values, but descendants can override these."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary20(self,
                          field = "Field<%(Dimension)s, std::vector<Scalar>>&"):
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary21(self,
                          field = "Field<%(Dimension)s, std::vector<Vector>>&"):
        return "void"

    @PYB11pycppname("applyGhostBoundary")
    @PYB11virtual
    @PYB11const
    def applyGhostBoundary22(self,
                             field = "Field<%(Dimension)s, RKCoefficients<%(Dimension)s>>&"):
        "We provide default copies for RKCoefficients, but descendants can override these."
        return "void"

    @PYB11pycppname("enforceBoundary")
    @PYB11virtual
    @PYB11const
    def enforceBoundary22(self,
                          field = "Field<%(Dimension)s, RKCoefficients<%(Dimension)s>>&"):
        "Default noop for RKCoefficients"
        return "void"

    @PYB11virtual
    def initializeProblemStartup(self,
                                 final = "const bool"):
        "Some boundaries need to know when a problem is starting up and all the physics packages have been initialized."
        return "void"

    @PYB11virtual
    @PYB11const
    def finalizeGhostBoundary(self):
        "Provide an optional hook that is to be called when all ghost boundaries are to have been set."
        return "void"

    @PYB11virtual
    def reset(self,
              dataBase = "const DataBase<%(Dimension)s>&"):
        "Overridable hook for clearing out the boundary condition."
        return "void"

    @PYB11virtual
    @PYB11const
    def clip(self,
             xmin = "Vector&",
             xmax = "Vector&"):
        "Optionally the boundary can clip an input box range. Defaults to no-op."
        return "void"

    @PYB11virtual
    @PYB11const
    def allowGhostCulling(self):
        "Optionally opt-out of ghost node culling."
        return "bool"

    @PYB11virtual
    @PYB11const
    def meshGhostNodes(self):
        "Some boundaries have ghosts we should exclude from tessellations. Provide a hook to note such cases."
        return "bool"

    @PYB11virtual
    @PYB11protected
    def addNodeList(self, nodeList="NodeList<%(Dimension)s>&"):
        return "void"

    #...........................................................................
    # Methods
    @PYB11const
    def haveNodeList(self, nodeList="const NodeList<%(Dimension)s>&"):
        "Check if we have entries for the given NodeList."
        return "bool"

    @PYB11returnpolicy("reference_internal")
    def accessBoundaryNodes(self, nodeList="NodeList<%(Dimension)s>&"):
        return "BoundaryNodes&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def controlNodes(self, nodeList="const NodeList<%(Dimension)s>&"):
        "Control nodes for a given NodeList"
        return "const std::vector<size_t>&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def ghostNodes(self, nodeList="const NodeList<%(Dimension)s>&"):
        "Ghost nodes for a given NodeList"
        return "const std::vector<size_t>&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def violationNodes(self, nodeList="const NodeList<%(Dimension)s>&"):
        "Violation nodes for a given NodeList"
        return "const std::vector<size_t>&"

    #...........................................................................
    # applyFieldListGhostBoundary/enforceFieldListBoundary
    @PYB11template("Value")
    @PYB11const
    def applyFieldListGhostBoundary(self, fieldList="FieldList<%(Dimension)s, %(Value)s>&"):
        "Apply ghost boundary to all Fields in FieldList"
        return "void"

    @PYB11template("Value")
    @PYB11const
    def enforceFieldListBoundary(self, fieldList="FieldList<%(Dimension)s, %(Value)s>&"):
        "Enforce boundary on all Fields in FieldList"
        return "void"

    for T, L in (("int", "int"),
                 ("uint64_t", "ULL"),
                 ("Scalar", "Scalar"),
                 ("Vector", "Vector"),
                 ("Tensor", "Tensor"),
                 ("SymTensor", "SymTensor"),
                 ("ThirdRankTensor", "ThirdRankTensor"),
                 ("FourthRankTensor", "FourthRankTensor"),
                 ("FifthRankTensor", "FifthRankTensor"),
                 ("std::vector<Scalar>", "StdVectorScalar"),
                 ("std::vector<Vector>", "StdVectorVector"),
                 ("RKCoefficients<%(Dimension)s>", "RKCoefficients")):
        exec('''
aflgb%(L)s = PYB11TemplateMethod(applyFieldListGhostBoundary, template_parameters="%(T)s", pyname="applyFieldListGhostBoundary")
eflgb%(L)s = PYB11TemplateMethod(enforceFieldListBoundary, template_parameters="%(T)s", pyname="enforceFieldListBoundary")
''' % {"T" : T,
       "L" : L})

    #...........................................................................
    # Attributes
    numGhostNodes = PYB11property("int")

    #---------------------------------------------------------------------------
    # BoundaryNodes
    #---------------------------------------------------------------------------
    class BoundaryNodes:

        controlNodes = PYB11readwrite()
        ghostNodes = PYB11readwrite()
        violationNodes = PYB11readwrite()

#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(BoundaryAbstractMethods, Boundary, pure_virtual=True)
