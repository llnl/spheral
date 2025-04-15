#-------------------------------------------------------------------------------
# InflowOutflowBoundary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Boundary import *
from Physics import *
from BoundaryAbstractMethods import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class InflowOutflowBoundary(Boundary, Physics):
    """InflowOutflowBoundary -- creates inflow ghost images, which become internal nodes
as they cross the specified boundary plane."""

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ThirdRankTensor = typename %(Dimension)s::ThirdRankTensor;
    using FourthRankTensor = typename %(Dimension)s::FourthRankTensor;
    using FifthRankTensor = typename %(Dimension)s::FifthRankTensor;
    using FacetedVolume = typename %(Dimension)s::FacetedVolume;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
    using Plane = GeomPlane<%(Dimension)s>;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               dataBase = "DataBase<%(Dimension)s>&",
               plane = "const GeomPlane<%(Dimension)s>&",
               empty = ("const bool", "false")):
        "Constructor"

    #...........................................................................
    # Methods
    @PYB11pycppname("applyGhostBoundary")
    @PYB11const
    def applyGhostBoundary0(self,
                            fieldBase = "FieldBase<%(Dimension)s>&"):
        "Apply the boundary condition to the ghost node values in the given Field."
        return "void"

    @PYB11virtual
    def cullGhostNodes(self,
                       flagSet = "const FieldList<%(Dimension)s, int>&",
                       old2newIndexMap = "FieldList<%(Dimension)s, int>&",
                       numNodesRemoved = "std::vector<int>&"):
        return "void"

    @PYB11virtual
    def initializeProblemStartup(self,
                                 final = "const bool"):
        "After physics have been initialized we take a snapshot of the node state."
        return "void"

    @PYB11virtual
    def finalize(self,
                 time = "const Scalar",
                 dt = "const Scalar",
                 dataBase = "DataBase<%(Dimension)s>&",
                 state = "State<%(Dimension)s>&",
                 derivs = "StateDerivatives<%(Dimension)s>&"):
        """Packages might want a hook to do some post-step finalizations.
Really we should rename this post-step finalize."""
        return "void"

    @PYB11const
    def numInflowNodes(self,
                       nodeList = "const NodeList<%(Dimension)s>&"):
        "Number of nodes in inflow stencil for the given NodeList"
        return "int"

    def clearStoredValues(self):
        "Clear out the stored values and ghost nodes"
        return "void"

    @PYB11virtual
    @PYB11const
    def allowGhostCulling(self):
        "Optionally opt-out of ghost node culling."
        return "bool"

    #...........................................................................
    @PYB11template("DataType")
    @PYB11returnpolicy("reference")
    def storedValues(self,
                     key = "const std::string",
                     dummy = ("const %(DataType)s&", "%(DataType)s()")):
        "Get the stored data for generating ghost nodes."
        return "std::vector<%(DataType)s>"

    storedValues_int =             PYB11TemplateMethod(storedValues, "int")
    storedValues_Scalar =          PYB11TemplateMethod(storedValues, "Scalar")
    storedValues_Vector =          PYB11TemplateMethod(storedValues, "Vector")
    storedValues_Tensor =          PYB11TemplateMethod(storedValues, "Tensor")
    storedValues_SymTensor =       PYB11TemplateMethod(storedValues, "SymTensor")
    storedValues_ThirdRankTensor = PYB11TemplateMethod(storedValues, "ThirdRankTensor")
    storedValues_FourthRankTensor = PYB11TemplateMethod(storedValues, "FourthRankTensor")
    storedValues_FifthRankTensor = PYB11TemplateMethod(storedValues, "FifthRankTensor")
    storedValues_FacetedVolume =   PYB11TemplateMethod(storedValues, "FacetedVolume")
    storedValues_vectorScalar =    PYB11TemplateMethod(storedValues, "std::vector<Scalar>")
    storedValues_vectorVector =    PYB11TemplateMethod(storedValues, "std::vector<Vector>")

    #...........................................................................
    @PYB11template("DataType")
    @PYB11returnpolicy("reference")
    @PYB11cppname("storedValues")
    def storedValuesF(self,
                      field = "const Field<%(Dimension)s, %(DataType)s>&"):
        "Get the stored data for generating ghost nodes."
        return "std::vector<%(DataType)s>"

    storedValuesF_int =             PYB11TemplateMethod(storedValuesF, "int", pyname="storedValues")
    storedValuesF_Scalar =          PYB11TemplateMethod(storedValuesF, "Scalar", pyname="storedValues")
    storedValuesF_Vector =          PYB11TemplateMethod(storedValuesF, "Vector", pyname="storedValues")
    storedValuesF_Tensor =          PYB11TemplateMethod(storedValuesF, "Tensor", pyname="storedValues")
    storedValuesF_SymTensor =       PYB11TemplateMethod(storedValuesF, "SymTensor", pyname="storedValues")
    storedValuesF_ThirdRankTensor = PYB11TemplateMethod(storedValuesF, "ThirdRankTensor", pyname="storedValues")
    storedValuesF_FourthRankTensor = PYB11TemplateMethod(storedValuesF, "FourthRankTensor", pyname="storedValues")
    storedValuesF_FifthRankTensor = PYB11TemplateMethod(storedValuesF, "FifthRankTensor", pyname="storedValues")
    storedValuesF_FacetedVolume =   PYB11TemplateMethod(storedValuesF, "FacetedVolume", pyname="storedValues")
    storedValuesF_vectorScalar =    PYB11TemplateMethod(storedValuesF, "std::vector<Scalar>", pyname="storedValues")
    storedValuesF_vectorVector =    PYB11TemplateMethod(storedValuesF, "std::vector<Vector>", pyname="storedValues")

    #...........................................................................
    @PYB11template("DataType")
    def setStoredValues(self,
                        key = "const std::string",
                        values = "const std::vector<%(DataType)s>&"):
        "Set the stored data for generating ghost nodes for a given Field key."
        return "void"

    setStoredValues_int =             PYB11TemplateMethod(setStoredValues, "int")
    setStoredValues_Scalar =          PYB11TemplateMethod(setStoredValues, "Scalar")
    setStoredValues_Vector =          PYB11TemplateMethod(setStoredValues, "Vector")
    setStoredValues_Tensor =          PYB11TemplateMethod(setStoredValues, "Tensor")
    setStoredValues_SymTensor =       PYB11TemplateMethod(setStoredValues, "SymTensor")
    setStoredValues_ThirdRankTensor = PYB11TemplateMethod(setStoredValues, "ThirdRankTensor")
    setStoredValues_FacetedVolume =   PYB11TemplateMethod(setStoredValues, "FacetedVolume")
    setStoredValues_vectorScalar =    PYB11TemplateMethod(setStoredValues, "std::vector<Scalar>")
    setStoredValues_vectorVector =    PYB11TemplateMethod(setStoredValues, "std::vector<Vector>")

    #...........................................................................
    @PYB11template("DataType")
    @PYB11cppname("setStoredValues")
    def setStoredValuesF(self,
                         field = "const Field<%(Dimension)s, %(DataType)s>&",
                         value = "const std::vector<%(DataType)s>&"):
        "Set the stored data for generating ghost nodes for a given field."
        return "void"

    setStoredValuesF_int =             PYB11TemplateMethod(setStoredValuesF, "int", pyname="setStoredValues")
    setStoredValuesF_Scalar =          PYB11TemplateMethod(setStoredValuesF, "Scalar", pyname="setStoredValues")
    setStoredValuesF_Vector =          PYB11TemplateMethod(setStoredValuesF, "Vector", pyname="setStoredValues")
    setStoredValuesF_Tensor =          PYB11TemplateMethod(setStoredValuesF, "Tensor", pyname="setStoredValues")
    setStoredValuesF_SymTensor =       PYB11TemplateMethod(setStoredValuesF, "SymTensor", pyname="setStoredValues")
    setStoredValuesF_ThirdRankTensor = PYB11TemplateMethod(setStoredValuesF, "ThirdRankTensor", pyname="setStoredValues")
    setStoredValuesF_FacetedVolume =   PYB11TemplateMethod(setStoredValuesF, "FacetedVolume", pyname="setStoredValues")
    setStoredValuesF_vectorScalar =    PYB11TemplateMethod(setStoredValuesF, "std::vector<Scalar>", pyname="setStoredValues")
    setStoredValuesF_vectorVector =    PYB11TemplateMethod(setStoredValuesF, "std::vector<Vector>", pyname="setStoredValues")

    #...........................................................................
    @PYB11template("DataType")
    @PYB11cppname("setStoredValues")
    def setConstantStoredValues(self,
                                key = "const std::string",
                                value = "const %(DataType)s&"):
        "Set the stored data for generating ghost nodes for a given Field key."
        return "void"

    setConstantStoredValues_int =             PYB11TemplateMethod(setConstantStoredValues, "int")
    setConstantStoredValues_Scalar =          PYB11TemplateMethod(setConstantStoredValues, "Scalar")
    setConstantStoredValues_Vector =          PYB11TemplateMethod(setConstantStoredValues, "Vector")
    setConstantStoredValues_Tensor =          PYB11TemplateMethod(setConstantStoredValues, "Tensor")
    setConstantStoredValues_SymTensor =       PYB11TemplateMethod(setConstantStoredValues, "SymTensor")
    setConstantStoredValues_ThirdRankTensor = PYB11TemplateMethod(setConstantStoredValues, "ThirdRankTensor")
    setConstantStoredValues_FacetedVolume =   PYB11TemplateMethod(setConstantStoredValues, "FacetedVolume")
    setConstantStoredValues_vectorScalar =    PYB11TemplateMethod(setConstantStoredValues, "std::vector<Scalar>")
    setConstantStoredValues_vectorVector =    PYB11TemplateMethod(setConstantStoredValues, "std::vector<Vector>")

    #...........................................................................
    @PYB11template("DataType")
    @PYB11cppname("setStoredValues")
    def setConstantStoredValuesF(self,
                                 field = "const Field<%(Dimension)s, %(DataType)s>&",
                                 value = "const %(DataType)s&"):
        "Set the stored data for generating ghost nodes for a given field."
        return "void"

    setConstantStoredValuesF_int =             PYB11TemplateMethod(setConstantStoredValuesF, "int", pyname="setConstantStoredValues")
    setConstantStoredValuesF_Scalar =          PYB11TemplateMethod(setConstantStoredValuesF, "Scalar", pyname="setConstantStoredValues")
    setConstantStoredValuesF_Vector =          PYB11TemplateMethod(setConstantStoredValuesF, "Vector", pyname="setConstantStoredValues")
    setConstantStoredValuesF_Tensor =          PYB11TemplateMethod(setConstantStoredValuesF, "Tensor", pyname="setConstantStoredValues")
    setConstantStoredValuesF_SymTensor =       PYB11TemplateMethod(setConstantStoredValuesF, "SymTensor", pyname="setConstantStoredValues")
    setConstantStoredValuesF_ThirdRankTensor = PYB11TemplateMethod(setConstantStoredValuesF, "ThirdRankTensor", pyname="setConstantStoredValues")
    setConstantStoredValuesF_FacetedVolume =   PYB11TemplateMethod(setConstantStoredValuesF, "FacetedVolume", pyname="setConstantStoredValues")
    setConstantStoredValuesF_vectorScalar =    PYB11TemplateMethod(setConstantStoredValuesF, "std::vector<Scalar>", pyname="setConstantStoredValues")
    setConstantStoredValuesF_vectorVector =    PYB11TemplateMethod(setConstantStoredValuesF, "std::vector<Vector>", pyname="setConstantStoredValues")

    #...........................................................................
    # Properties
    inflowRadius = PYB11property("Scalar", "inflowRadius", "inflowRadius",doc="radius of inflow clip")
    dataBase = PYB11property(doc="The DataBase for the NodeLists we know about")
    plane = PYB11property(doc="The inflowplane")
    storedKeys = PYB11property(doc="Keys for all the Fields we have stored ghost information about")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, InflowOutflowBoundary)
PYB11inject(BoundaryAbstractMethods, InflowOutflowBoundary, virtual=True, pure_virtual=False)
PYB11inject(PhysicsAbstractMethods, InflowOutflowBoundary, virtual=True, pure_virtual=False)
