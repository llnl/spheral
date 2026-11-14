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
  using Scalar = %(Dimension)s::Scalar;
  using Vector = %(Dimension)s::Vector;
  using Tensor = %(Dimension)s::Tensor;
  using SymTensor = %(Dimension)s::SymTensor;
  using ThirdRankTensor = %(Dimension)s::ThirdRankTensor;
  using FourthRankTensor = %(Dimension)s::FourthRankTensor;
  using FifthRankTensor = %(Dimension)s::FifthRankTensor;
  using FacetedVolume = %(Dimension)s::FacetedVolume;
  using Plane = GeomPlane<%(Dimension)s>;
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
