#-------------------------------------------------------------------------------
# SphericalOriginBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ReflectingBoundary import *
from BoundaryAbstractMethods import *

@PYB11template()
@PYB11template_dict({"Dimension" : "Dim<1>"})
class SphericalOriginBoundary(ReflectingBoundary):
    PYB11typedefs = """
    typedef %(Dimension)s::Scalar Scalar;
    typedef %(Dimension)s::Vector Vector;
    typedef %(Dimension)s::Tensor Tensor;
    typedef %(Dimension)s::SymTensor SymTensor;
    typedef %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef %(Dimension)s::FourthRankTensor FourthRankTensor;
    typedef %(Dimension)s::FifthRankTensor FifthRankTensor;
    typedef %(Dimension)s::FacetedVolume FacetedVolume;
"""

    #...........................................................................
    # Constructors
    def pyinit(self, dataBase="const DataBase<%(Dimension)s>&"):
        "Construct with the DataBase"

    #...........................................................................
    # Methods
    @PYB11virtual
    def setGhostNodes(self,
                      nodeList = "NodeList<%(Dimension)s>&"):
        return "void"

    @PYB11virtual
    def updateGhostNodes(self,
                         nodeList = "NodeList<%(Dimension)s>&"):
        return "void"


    @PYB11virtual
    def setViolationNodes(self, nodeList="NodeList<%(Dimension)s>&"):
        return "void"

    @PYB11virtual
    def updateViolationNodes(self, nodeList="NodeList<%(Dimension)s>&"):
        return "void"

    @PYB11virtual
    @PYB11const
    def label(self):
        "The label for writing in restart files"
        return "std::string"

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
#PYB11inject(BoundaryAbstractMethods, SphericalOriginBoundary, virtual=True, pure_virtual=False)
