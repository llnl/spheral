#-------------------------------------------------------------------------------
# DataBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from DataBase import *

@PYB11template("Dimension")
@PYB11dynamic_attr
class DataBase:

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename %(Dimension)s::FourthRankTensor FourthRankTensor;
    typedef typename %(Dimension)s::FifthRankTensor FifthRankTensor;
    typedef typename %(Dimension)s::FacetedVolume FacetedVolume;
    typedef typename DataBase<%(Dimension)s>::ConnectivityMapPtr ConnectivityMapPtr;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    #...........................................................................
    # Methods
    @PYB11const
    def reinitializeNeighbors(self):
        "Optimize all Neighbor objects for the current state."
        return "void"

    @PYB11const
    def updateConnectivityMap(self,
                              computeGhostConnectivity = ("const bool", "false"),
                              computeOverlapConnectivity = ("const bool", "false"),
                              computeIntersectionConnectivity = ("const bool", "false")):
        "Update the internal connectivity map."
        return "void"

    @PYB11const
    def patchConnectivityMap(self,
                             flags = "const FieldList<%(Dimension)s, size_t>&",
                             old2new = "const FieldList<%(Dimension)s, size_t>&"):
        "Update the internal connectivity map."
        return "void"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def connectivityMap(self,
                        computeGhostConnectivity=("const bool", "false"),
                        computeOverlapConnectivity=("const bool", "false"),
                        computeIntersectionConnectivity=("const bool", "false")):
        "Get the connectivity map, optionally including ghost, overlap, or intersection connectivity"
        return "const ConnectivityMap<%(Dimension)s>&"

    @PYB11const
    def connectivityMapPtr(self,
                           computeGhostConnectivity=("const bool", "false"),
                           computeOverlapConnectivity=("const bool", "false"),
                           computeIntersectionConnectivity=("const bool", "false")):
        "Get the connectivity map as a std::shared_ptr, optionally including ghost, overlap, or intersection connectivity"
        return "ConnectivityMapPtr"

    def appendNodeList(self, nodeList="SolidNodeList<%(Dimension)s>&"):
        "Add a SolidNodeList"
        return "void"

    @PYB11pycppname("appendNodeList")
    def appendNodeList1(self, nodeList="FluidNodeList<%(Dimension)s>&"):
        "Add a FluidNodeList"
        return "void"

    @PYB11pycppname("appendNodeList")
    def appendNodeList2(self, nodeList="NodeList<%(Dimension)s>&"):
        "Add a NodeList"
        return "void"

    @PYB11pycppname("appendNodeList")
    def appendNodeList2(self, nodeList="DEMNodeList<%(Dimension)s>&"):
        "Add a DEMNodeList"
        return "void"

    def deleteNodeList(self, nodeList="SolidNodeList<%(Dimension)s>&"):
        "Remove a SolidNodeList"
        return "void"

    @PYB11pycppname("deleteNodeList")
    def deleteNodeList1(self, nodeList="FluidNodeList<%(Dimension)s>&"):
        "Remove a FluidNodeList"
        return "void"

    @PYB11pycppname("deleteNodeList")
    def deleteNodeList2(self, nodeList="NodeList<%(Dimension)s>&"):
        "Remove a NodeList"
        return "void"

    @PYB11pycppname("deleteNodeList")
    def deleteNodeList2(self, nodeList="DEMNodeList<%(Dimension)s>&"):
        "Remove a DEMNodeList"
        return "void"

    @PYB11const
    def haveNodeList(self, nodeList="const NodeList<%(Dimension)s>&"):
        "Check if a NodeList is in the DataBase"
        return "bool"

    @PYB11const
    def setMasterNodeLists(self,
                           position = "const Vector&",
                           H = "const SymTensor&",
                           masterLists = "std::vector<std::vector<int>>&",
                           coarseNeighbors = "std::vector<std::vector<int>>&",
                           computeGhostConnectivity = "bool"):
        "Set the master/coarse neighbor lists for all NodeLists"
        return "void"

    @PYB11const
    def setMasterFluidNodeLists(self,
                                position = "const Vector&",
                                H = "const SymTensor&",
                                masterLists = "std::vector<std::vector<int>>&",
                                coarseNeighbors = "std::vector<std::vector<int>>&",
                                computeGhostConnectivity = "bool"):
        "Set the master/coarse neighbor lists for all NodeLists"
        return "void"

    @PYB11const
    def setRefineNodeLists(self,
                           position = "const Vector&",
                           H = "const SymTensor&",
                           coarseNeighbors = "const std::vector<std::vector<int>>&",
                           refineNeighbors = "std::vector<std::vector<int>>&"):
        "Set the refine neighbor lists for all NodeLists"
        return "void"

    @PYB11const
    def setRefineFluidNodeLists(self,
                                position = "const Vector&",
                                H = "const SymTensor&",
                                coarseNeighbors = "const std::vector<std::vector<int>>&",
                                refineNeighbors = "std::vector<std::vector<int>>&"):
        "Set the refine neighbor lists for all FluidNodeLists"
        return "void"

    @PYB11const
    def globalHinverse(self, result="FieldList<%(Dimension)s, SymTensor>&"):
        return "void"

    @PYB11const
    def fluidHinverse(self, result="FieldList<%(Dimension)s, SymTensor>&"):
        return "void"

    @PYB11const
    def fluidPressure(self, result="FieldList<%(Dimension)s, Scalar>&"):
        return "void"

    @PYB11const
    def fluidTemperature(self, result="FieldList<%(Dimension)s, Scalar>&"):
        return "void"

    @PYB11const
    def fluidSoundSpeed(self, result="FieldList<%(Dimension)s, Scalar>&"):
        return "void"

    @PYB11const
    def fluidVolume(self, result="FieldList<%(Dimension)s, Scalar>&"):
        return "void"

    @PYB11const
    def fluidGamma(self, result="FieldList<%(Dimension)s, Scalar>&"):
        return "void"

    @PYB11const
    def fluidEntropy(self, result="FieldList<%(Dimension)s, Scalar>&"):
        return "void"

    @PYB11const
    def fluidLinearMomentum(self, result="FieldList<%(Dimension)s, Vector>&"):
        return "void"

    @PYB11const
    def fluidTotalEnergy(self, result="FieldList<%(Dimension)s, Scalar>&"):
        return "void"

    @PYB11const
    def fluidSpecificHeat(self,
                          temperature = "const FieldList<%(Dimension)s, Scalar>&",
                          result="FieldList<%(Dimension)s, Scalar>&"):
        return "void"
    
    @PYB11const
    def boundingBox(self,
                    xmin = "Vector&",
                    xmax = "Vector&",
                    ghost = ("const bool", "true")):
        "Compute coordinates bounding all nodes in the DataBase."
        return "void"

    @PYB11const
    @PYB11pycppname("boundingBox")
    def boundingBox1(self,
                     xmin = "Vector&",
                     xmax = "Vector&",
                     mask = "const FieldList<%(Dimension)s, int>&",
                     ghost = ("const bool", "true")):
        "Compute coordinates bounding all nodes in the DataBase."
        return "void"

    @PYB11const
    def localSamplingBoundingVolume(self,
                                    centroid = "Vector&",
                                    radiusNodes = "double&",
                                    radiusSample = "double&",
                                    xminNodes = "Vector&",
                                    xmaxNodes = "Vector&",
                                    xminSample = "Vector&",
                                    xmaxSample = "Vector&"):
        "Return the local max sampling extents."
        return "void"

    @PYB11const
    def globalSamplingBoundingVolume(self,
                                     centroid = "Vector&",
                                     radiusNodes = "double&",
                                     radiusSample = "double&",
                                     xminNodes = "Vector&",
                                     xmaxNodes = "Vector&",
                                     xminSample = "Vector&",
                                     xmaxSample = "Vector&"):
        "Return the global max sampling extents."
        return "void"

    @PYB11const
    def localSamplingBoundingBoxes(self,
                                   xminima = "std::vector<Vector>&",
                                   xmaxima = "std::vector<Vector>&"):
        "Return the local min and max sampling extents for groupings of connected nodes."
        return "void"

    @PYB11const
    def globalSamplingBoundingBoxes(self,
                                    xminima = "std::vector<Vector>&",
                                    xmaxima = "std::vector<Vector>&"):
        "Return the global min and max sampling extents for groupings of connected nodes."
        return "void"

    @PYB11const
    def valid(self):
        "Provide a method to determine if the DataBase is in a minimally defined valid state."
        return "bool"

    #...........................................................................
    # FieldList generation methods
    @PYB11template("DataType")
    @PYB11const
    def newGlobalFieldList(self,
                           value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()"),
                           name = ("const Field<%(Dimension)s, %(DataType)s>::FieldName", '"Unnamed Field"')):
        "Construct a new FieldList<%(DataType)s> for all NodeLists in DataBase"
        return "FieldList<%(Dimension)s, %(DataType)s>"

    @PYB11template("DataType")
    @PYB11const
    def newFluidFieldList(self,
                          value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()"),
                          name = ("const Field<%(Dimension)s, %(DataType)s>::FieldName", '"Unnamed Field"')):
        "Construct a new FieldList<%(DataType)s> for all FluidNodeLists in DataBase"
        return "FieldList<%(Dimension)s, %(DataType)s>"

    @PYB11template("DataType")
    @PYB11const
    def newSolidFieldList(self,
                          value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()"),
                          name = ("const Field<%(Dimension)s, %(DataType)s>::FieldName", '"Unnamed Field"')):
        "Construct a new FieldList<%(DataType)s> for all SolidNodeLists in DataBase"
        return "FieldList<%(Dimension)s, %(DataType)s>"

    @PYB11template("DataType")
    @PYB11const
    def newDEMFieldList(self,
                        value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()"),
                        name = ("const Field<%(Dimension)s, %(DataType)s>::FieldName", '"Unnamed Field"')):
        "Construct a new FieldList<%(DataType)s> for all DEMNodeLists in DataBase"
        return "FieldList<%(Dimension)s, %(DataType)s>"

    @PYB11template("DataType")
    @PYB11const
    def resizeGlobalFieldList(self,
                              fieldList = "FieldList<%(Dimension)s, %(DataType)s>&",
                              value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()"),
                              name = ("const Field<%(Dimension)s, %(DataType)s>::FieldName", '"Unnamed Field"'),
                              resetValues = ("const bool", "true")):
        """Resize a FieldList to the number of NodeLists.
Optionally we can also set all elements in the FieldList to the specified value.
Note that if the FieldList is resized it is reconstructed from scratch, so all elements
will get the new value regardless of resetValues."""
        return "void"

    @PYB11template("DataType")
    @PYB11const
    def resizeFluidFieldList(self,
                             fieldList = "FieldList<%(Dimension)s, %(DataType)s>&",
                             value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()"),
                             name = ("const Field<%(Dimension)s, %(DataType)s>::FieldName", '"Unnamed Field"'),
                             resetValues = ("const bool", "true")):
        """Resize a FieldList to the number of FluidNodeLists.
Optionally we can also set all elements in the FieldList to the specified value.
Note that if the FieldList is resized it is reconstructed from scratch, so all elements
will get the new value regardless of resetValues."""
        return "void"

    @PYB11template("DataType")
    @PYB11const
    def resizeSolidFieldList(self,
                             fieldList = "FieldList<%(Dimension)s, %(DataType)s>&",
                             value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()"),
                             name = ("const Field<%(Dimension)s, %(DataType)s>::FieldName", '"Unnamed Field"'),
                             resetValues = ("const bool", "true")):
        """Resize a FieldList to the number of SolidNodeLists.
Optionally we can also set all elements in the FieldList to the specified value.
Note that if the FieldList is resized it is reconstructed from scratch, so all elements
will get the new value regardless of resetValues."""
        return "void"

    @PYB11template("DataType")
    @PYB11const
    def resizeDEMFieldList(self,
                             fieldList = "FieldList<%(Dimension)s, %(DataType)s>&",
                             value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()"),
                             name = ("const Field<%(Dimension)s, %(DataType)s>::FieldName", '"Unnamed Field"'),
                             resetValues = ("const bool", "true")):
        """Resize a FieldList to the number of DEMNodeLists.
Optionally we can also set all elements in the FieldList to the specified value.
Note that if the FieldList is resized it is reconstructed from scratch, so all elements
will get the new value regardless of resetValues."""
        return "void"

    newGlobalIntFieldList              = PYB11TemplateMethod(newGlobalFieldList, template_parameters="int")
    newGlobalScalarFieldList           = PYB11TemplateMethod(newGlobalFieldList, template_parameters="double")
    newGlobalVectorFieldList           = PYB11TemplateMethod(newGlobalFieldList, template_parameters="Vector")
    newGlobalTensorFieldList           = PYB11TemplateMethod(newGlobalFieldList, template_parameters="Tensor")
    newGlobalSymTensorFieldList        = PYB11TemplateMethod(newGlobalFieldList, template_parameters="SymTensor")
    newGlobalThirdRankTensorFieldList  = PYB11TemplateMethod(newGlobalFieldList, template_parameters="ThirdRankTensor")
    newGlobalFourthRankTensorFieldList = PYB11TemplateMethod(newGlobalFieldList, template_parameters="FourthRankTensor")
    newGlobalFifthRankTensorFieldList  = PYB11TemplateMethod(newGlobalFieldList, template_parameters="FifthRankTensor")
    newGlobalFacetedVolumeFieldList    = PYB11TemplateMethod(newGlobalFieldList, template_parameters="FacetedVolume")
    newGlobalvector_of_intFieldList    = PYB11TemplateMethod(newGlobalFieldList, template_parameters="std::vector<int>")
    newGlobalvector_of_doubleFieldList = PYB11TemplateMethod(newGlobalFieldList, template_parameters="std::vector<double>")
    newGlobalvector_of_VectorFieldList = PYB11TemplateMethod(newGlobalFieldList, template_parameters="std::vector<Vector>")
    newGlobalvector_of_CellFaceFlagFieldList = PYB11TemplateMethod(newGlobalFieldList, template_parameters="std::vector<CellFaceFlag>")
    newGlobalDomainNodeFieldList       = PYB11TemplateMethod(newGlobalFieldList, template_parameters="DomainNode<%(Dimension)s>")
    
    newFluidIntFieldList              = PYB11TemplateMethod(newFluidFieldList, template_parameters="int")
    newFluidScalarFieldList           = PYB11TemplateMethod(newFluidFieldList, template_parameters="double")
    newFluidVectorFieldList           = PYB11TemplateMethod(newFluidFieldList, template_parameters="Vector")
    newFluidTensorFieldList           = PYB11TemplateMethod(newFluidFieldList, template_parameters="Tensor")
    newFluidSymTensorFieldList        = PYB11TemplateMethod(newFluidFieldList, template_parameters="SymTensor")
    newFluidThirdRankTensorFieldList  = PYB11TemplateMethod(newFluidFieldList, template_parameters="ThirdRankTensor")
    newFluidFourthRankTensorFieldList = PYB11TemplateMethod(newFluidFieldList, template_parameters="FourthRankTensor")
    newFluidFifthRankTensorFieldList  = PYB11TemplateMethod(newFluidFieldList, template_parameters="FifthRankTensor")
    newFluidFacetedVolumeFieldList    = PYB11TemplateMethod(newFluidFieldList, template_parameters="FacetedVolume")
    newFluidvector_of_intFieldList    = PYB11TemplateMethod(newFluidFieldList, template_parameters="std::vector<int>")
    newFluidvector_of_doubleFieldList = PYB11TemplateMethod(newFluidFieldList, template_parameters="std::vector<double>")
    newFluidvector_of_VectorFieldList = PYB11TemplateMethod(newFluidFieldList, template_parameters="std::vector<Vector>")
    newFluidvector_of_CellFaceFlagFieldList = PYB11TemplateMethod(newFluidFieldList, template_parameters="std::vector<CellFaceFlag>")
    newFluidDomainNodeFieldList       = PYB11TemplateMethod(newFluidFieldList, template_parameters="DomainNode<%(Dimension)s>")

    newSolidIntFieldList              = PYB11TemplateMethod(newSolidFieldList, template_parameters="int")
    newSolidScalarFieldList           = PYB11TemplateMethod(newSolidFieldList, template_parameters="double")
    newSolidVectorFieldList           = PYB11TemplateMethod(newSolidFieldList, template_parameters="Vector")
    newSolidTensorFieldList           = PYB11TemplateMethod(newSolidFieldList, template_parameters="Tensor")
    newSolidSymTensorFieldList        = PYB11TemplateMethod(newSolidFieldList, template_parameters="SymTensor")
    newSolidThirdRankTensorFieldList  = PYB11TemplateMethod(newSolidFieldList, template_parameters="ThirdRankTensor")
    newSolidFourthRankTensorFieldList = PYB11TemplateMethod(newSolidFieldList, template_parameters="FourthRankTensor")
    newSolidFifthRankTensorFieldList  = PYB11TemplateMethod(newSolidFieldList, template_parameters="FifthRankTensor")
    newSolidFacetedVolumeFieldList    = PYB11TemplateMethod(newSolidFieldList, template_parameters="FacetedVolume")
    newSolidvector_of_intFieldList    = PYB11TemplateMethod(newSolidFieldList, template_parameters="std::vector<int>")
    newSolidvector_of_doubleFieldList = PYB11TemplateMethod(newSolidFieldList, template_parameters="std::vector<double>")
    newSolidvector_of_VectorFieldList = PYB11TemplateMethod(newSolidFieldList, template_parameters="std::vector<Vector>")
    newSolidvector_of_CellFaceFlagFieldList = PYB11TemplateMethod(newSolidFieldList, template_parameters="std::vector<CellFaceFlag>")
    newSolidDomainNodeFieldList       = PYB11TemplateMethod(newSolidFieldList, template_parameters="DomainNode<%(Dimension)s>")

    newDEMIntFieldList              = PYB11TemplateMethod(newDEMFieldList, template_parameters="int")
    newDEMScalarFieldList           = PYB11TemplateMethod(newDEMFieldList, template_parameters="double")
    newDEMVectorFieldList           = PYB11TemplateMethod(newDEMFieldList, template_parameters="Vector")
    newDEMTensorFieldList           = PYB11TemplateMethod(newDEMFieldList, template_parameters="Tensor")
    newDEMSymTensorFieldList        = PYB11TemplateMethod(newDEMFieldList, template_parameters="SymTensor")
    newDEMThirdRankTensorFieldList  = PYB11TemplateMethod(newDEMFieldList, template_parameters="ThirdRankTensor")
    newDEMFourthRankTensorFieldList = PYB11TemplateMethod(newDEMFieldList, template_parameters="FourthRankTensor")
    newDEMFifthRankTensorFieldList  = PYB11TemplateMethod(newDEMFieldList, template_parameters="FifthRankTensor")
    newDEMFacetedVolumeFieldList    = PYB11TemplateMethod(newDEMFieldList, template_parameters="FacetedVolume")
    newDEMvector_of_intFieldList    = PYB11TemplateMethod(newDEMFieldList, template_parameters="std::vector<int>")
    newDEMvector_of_doubleFieldList = PYB11TemplateMethod(newDEMFieldList, template_parameters="std::vector<double>")
    newDEMvector_of_VectorFieldList = PYB11TemplateMethod(newDEMFieldList, template_parameters="std::vector<Vector>")
    newDEMvector_of_CellFaceFlagFieldList = PYB11TemplateMethod(newDEMFieldList, template_parameters="std::vector<CellFaceFlag>")
    newDEMvector_of_CellFaceFlagFieldList = PYB11TemplateMethod(newDEMFieldList, template_parameters="std::vector<int>")

    resizeGlobalIntFieldList              = PYB11TemplateMethod(resizeGlobalFieldList, template_parameters="int")
    resizeGlobalScalarFieldList           = PYB11TemplateMethod(resizeGlobalFieldList, template_parameters="double")
    resizeGlobalVectorFieldList           = PYB11TemplateMethod(resizeGlobalFieldList, template_parameters="Vector")
    resizeGlobalTensorFieldList           = PYB11TemplateMethod(resizeGlobalFieldList, template_parameters="Tensor")
    resizeGlobalSymTensorFieldList        = PYB11TemplateMethod(resizeGlobalFieldList, template_parameters="SymTensor")
    resizeGlobalThirdRankTensorFieldList  = PYB11TemplateMethod(resizeGlobalFieldList, template_parameters="ThirdRankTensor")
    resizeGlobalFourthRankTensorFieldList = PYB11TemplateMethod(resizeGlobalFieldList, template_parameters="FourthRankTensor")
    resizeGlobalFifthRankTensorFieldList  = PYB11TemplateMethod(resizeGlobalFieldList, template_parameters="FifthRankTensor")
    resizeGlobalFacetedVolumeFieldList    = PYB11TemplateMethod(resizeGlobalFieldList, template_parameters="FacetedVolume")
    resizeGlobalvector_of_intFieldList    = PYB11TemplateMethod(resizeGlobalFieldList, template_parameters="std::vector<int>")
    resizeGlobalvector_of_doubleFieldList = PYB11TemplateMethod(resizeGlobalFieldList, template_parameters="std::vector<double>")
    resizeGlobalvector_of_VectorFieldList = PYB11TemplateMethod(resizeGlobalFieldList, template_parameters="std::vector<Vector>")
    resizeGlobalDomainNodeFieldList       = PYB11TemplateMethod(resizeGlobalFieldList, template_parameters="DomainNode<%(Dimension)s>")

    resizeFluidIntFieldList              = PYB11TemplateMethod(resizeFluidFieldList, template_parameters="int")
    resizeFluidScalarFieldList           = PYB11TemplateMethod(resizeFluidFieldList, template_parameters="double")
    resizeFluidVectorFieldList           = PYB11TemplateMethod(resizeFluidFieldList, template_parameters="Vector")
    resizeFluidTensorFieldList           = PYB11TemplateMethod(resizeFluidFieldList, template_parameters="Tensor")
    resizeFluidSymTensorFieldList        = PYB11TemplateMethod(resizeFluidFieldList, template_parameters="SymTensor")
    resizeFluidThirdRankTensorFieldList  = PYB11TemplateMethod(resizeFluidFieldList, template_parameters="ThirdRankTensor")
    resizeFluidFourthRankTensorFieldList = PYB11TemplateMethod(resizeFluidFieldList, template_parameters="FourthRankTensor")
    resizeFluidFifthRankTensorFieldList  = PYB11TemplateMethod(resizeFluidFieldList, template_parameters="FifthRankTensor")
    resizeFluidFacetedVolumeFieldList    = PYB11TemplateMethod(resizeFluidFieldList, template_parameters="FacetedVolume")
    resizeFluidvector_of_intFieldList    = PYB11TemplateMethod(resizeFluidFieldList, template_parameters="std::vector<int>")
    resizeFluidvector_of_doubleFieldList = PYB11TemplateMethod(resizeFluidFieldList, template_parameters="std::vector<double>")
    resizeFluidvector_of_VectorFieldList = PYB11TemplateMethod(resizeFluidFieldList, template_parameters="std::vector<Vector>")
    resizeFluidDomainNodeFieldList       = PYB11TemplateMethod(resizeFluidFieldList, template_parameters="DomainNode<%(Dimension)s>")

    resizeSolidIntFieldList              = PYB11TemplateMethod(resizeSolidFieldList, template_parameters="int")
    resizeSolidScalarFieldList           = PYB11TemplateMethod(resizeSolidFieldList, template_parameters="double")
    resizeSolidVectorFieldList           = PYB11TemplateMethod(resizeSolidFieldList, template_parameters="Vector")
    resizeSolidTensorFieldList           = PYB11TemplateMethod(resizeSolidFieldList, template_parameters="Tensor")
    resizeSolidSymTensorFieldList        = PYB11TemplateMethod(resizeSolidFieldList, template_parameters="SymTensor")
    resizeSolidThirdRankTensorFieldList  = PYB11TemplateMethod(resizeSolidFieldList, template_parameters="ThirdRankTensor")
    resizeSolidFourthRankTensorFieldList = PYB11TemplateMethod(resizeSolidFieldList, template_parameters="FourthRankTensor")
    resizeSolidFifthRankTensorFieldList  = PYB11TemplateMethod(resizeSolidFieldList, template_parameters="FifthRankTensor")
    resizeSolidFacetedVolumeFieldList    = PYB11TemplateMethod(resizeSolidFieldList, template_parameters="FacetedVolume")
    resizeSolidvector_of_intFieldList    = PYB11TemplateMethod(resizeSolidFieldList, template_parameters="std::vector<int>")
    resizeSolidvector_of_doubleFieldList = PYB11TemplateMethod(resizeSolidFieldList, template_parameters="std::vector<double>")
    resizeSolidDomainNodeFieldList       = PYB11TemplateMethod(resizeSolidFieldList, template_parameters="DomainNode<%(Dimension)s>")

    resizeDEMIntFieldList              = PYB11TemplateMethod(resizeDEMFieldList, template_parameters="int")
    resizeDEMScalarFieldList           = PYB11TemplateMethod(resizeDEMFieldList, template_parameters="double")
    resizeDEMVectorFieldList           = PYB11TemplateMethod(resizeDEMFieldList, template_parameters="Vector")
    resizeDEMTensorFieldList           = PYB11TemplateMethod(resizeDEMFieldList, template_parameters="Tensor")
    resizeDEMSymTensorFieldList        = PYB11TemplateMethod(resizeDEMFieldList, template_parameters="SymTensor")
    resizeDEMThirdRankTensorFieldList  = PYB11TemplateMethod(resizeDEMFieldList, template_parameters="ThirdRankTensor")
    resizeDEMFourthRankTensorFieldList = PYB11TemplateMethod(resizeDEMFieldList, template_parameters="FourthRankTensor")
    resizeDEMFifthRankTensorFieldList  = PYB11TemplateMethod(resizeDEMFieldList, template_parameters="FifthRankTensor")
    resizeDEMFacetedVolumeFieldList    = PYB11TemplateMethod(resizeDEMFieldList, template_parameters="FacetedVolume")
    resizeDEMvector_of_intFieldList    = PYB11TemplateMethod(resizeDEMFieldList, template_parameters="std::vector<int>")
    resizeDEMvector_of_doubleFieldList = PYB11TemplateMethod(resizeDEMFieldList, template_parameters="std::vector<double>")
    resizeDEMvector_of_VectorFieldList = PYB11TemplateMethod(resizeDEMFieldList, template_parameters="std::vector<Vector>")
    #...........................................................................
    # Array generation methods
    @PYB11template("DataType")
    @PYB11const
    def newGlobalArray(self,
                       value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()")):
        "Construct a new array<%(DataType)s> for all NodeLists in DataBase"
        return "std::vector<std::vector<%(DataType)s>>"

    @PYB11template("DataType")
    @PYB11const
    def newFluidArray(self,
                      value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()")):
        "Construct a new array<%(DataType)s> for all FluidNodeLists in DataBase"
        return "std::vector<std::vector<%(DataType)s>>"

    @PYB11template("DataType")
    @PYB11const
    def newSolidArray(self,
                      value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()")):
        "Construct a new array<%(DataType)s> for all SolidNodeLists in DataBase"
        return "std::vector<std::vector<%(DataType)s>>"

    @PYB11template("DataType")
    @PYB11const
    def newDEMArray(self,
                      value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()")):
        "Construct a new array<%(DataType)s> for all DEMNodeLists in DataBase"
        return "std::vector<std::vector<%(DataType)s>>"

    @PYB11template("DataType")
    @PYB11const
    def resizeGlobalArray(self,
                          array = "std::vector<std::vector<%(DataType)s>>&",
                          value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()"),
                          resetValues = ("const bool", "true")):
        """Resize an array to the number of NodeLists.
Optionally we can also set all elements in the Array to the specified value.
Note that if the Array is resized it is reconstructed from scratch, so all elements
will get the new value regardless of resetValues."""
        return "void"

    @PYB11template("DataType")
    @PYB11const
    def resizeFluidArray(self,
                         array = "std::vector<std::vector<%(DataType)s>>&",
                         value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()"),
                         resetValues = ("const bool", "true")):
        """Resize an array to the number of FluidNodeLists.
Optionally we can also set all elements in the Array to the specified value.
Note that if the Array is resized it is reconstructed from scratch, so all elements
will get the new value regardless of resetValues."""
        return "void"

    @PYB11template("DataType")
    @PYB11const
    def resizeSolidArray(self,
                         array = "std::vector<std::vector<%(DataType)s>>&",
                         value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()"),
                         resetValues = ("const bool", "true")):
        """Resize an array to the number of SolidNodeLists.
Optionally we can also set all elements in the Array to the specified value.
Note that if the Array is resized it is reconstructed from scratch, so all elements
will get the new value regardless of resetValues."""
        return "void"

    @PYB11template("DataType")
    @PYB11const
    def resizeDEMArray(self,
                         array = "std::vector<std::vector<%(DataType)s>>&",
                         value = ("const %(DataType)s", "DataTypeTraits<%(DataType)s>::zero()"),
                         resetValues = ("const bool", "true")):
        """Resize an array to the number of DEMNodeLists.
Optionally we can also set all elements in the Array to the specified value.
Note that if the Array is resized it is reconstructed from scratch, so all elements
will get the new value regardless of resetValues."""
        return "void"

    newGlobalIntArray              = PYB11TemplateMethod(newGlobalArray, template_parameters="int")
    newGlobalScalarArray           = PYB11TemplateMethod(newGlobalArray, template_parameters="double")
    newGlobalVectorArray           = PYB11TemplateMethod(newGlobalArray, template_parameters="Vector")
    newGlobalTensorArray           = PYB11TemplateMethod(newGlobalArray, template_parameters="Tensor")
    newGlobalSymTensorArray        = PYB11TemplateMethod(newGlobalArray, template_parameters="SymTensor")
    newGlobalThirdRankTensorArray  = PYB11TemplateMethod(newGlobalArray, template_parameters="ThirdRankTensor")
    newGlobalFourthRankTensorArray = PYB11TemplateMethod(newGlobalArray, template_parameters="FourthRankTensor")
    newGlobalFifthRankTensorArray  = PYB11TemplateMethod(newGlobalArray, template_parameters="FifthRankTensor")
    newGlobalFacetedVolumeArray    = PYB11TemplateMethod(newGlobalArray, template_parameters="FacetedVolume")
    newGlobalvector_of_intArray    = PYB11TemplateMethod(newGlobalArray, template_parameters="std::vector<int>")
    newGlobalvector_of_doubleArray = PYB11TemplateMethod(newGlobalArray, template_parameters="std::vector<double>")
    newGlobalvector_of_VectorArray = PYB11TemplateMethod(newGlobalArray, template_parameters="std::vector<Vector>")
    newGlobalvector_of_CellFaceFlagArray = PYB11TemplateMethod(newGlobalArray, template_parameters="std::vector<CellFaceFlag>")
    newGlobalDomainNodeArray       = PYB11TemplateMethod(newGlobalArray, template_parameters="DomainNode<%(Dimension)s>")
    
    newFluidIntArray              = PYB11TemplateMethod(newFluidArray, template_parameters="int")
    newFluidScalarArray           = PYB11TemplateMethod(newFluidArray, template_parameters="double")
    newFluidVectorArray           = PYB11TemplateMethod(newFluidArray, template_parameters="Vector")
    newFluidTensorArray           = PYB11TemplateMethod(newFluidArray, template_parameters="Tensor")
    newFluidSymTensorArray        = PYB11TemplateMethod(newFluidArray, template_parameters="SymTensor")
    newFluidThirdRankTensorArray  = PYB11TemplateMethod(newFluidArray, template_parameters="ThirdRankTensor")
    newFluidFourthRankTensorArray = PYB11TemplateMethod(newFluidArray, template_parameters="FourthRankTensor")
    newFluidFifthRankTensorArray  = PYB11TemplateMethod(newFluidArray, template_parameters="FifthRankTensor")
    newFluidFacetedVolumeArray    = PYB11TemplateMethod(newFluidArray, template_parameters="FacetedVolume")
    newFluidvector_of_intArray    = PYB11TemplateMethod(newFluidArray, template_parameters="std::vector<int>")
    newFluidvector_of_doubleArray = PYB11TemplateMethod(newFluidArray, template_parameters="std::vector<double>")
    newFluidvector_of_VectorArray = PYB11TemplateMethod(newFluidArray, template_parameters="std::vector<Vector>")
    newFluidvector_of_CellFaceFlagArray = PYB11TemplateMethod(newFluidArray, template_parameters="std::vector<CellFaceFlag>")
    newFluidDomainNodeArray       = PYB11TemplateMethod(newFluidArray, template_parameters="DomainNode<%(Dimension)s>")

    newSolidIntArray              = PYB11TemplateMethod(newSolidArray, template_parameters="int")
    newSolidScalarArray           = PYB11TemplateMethod(newSolidArray, template_parameters="double")
    newSolidVectorArray           = PYB11TemplateMethod(newSolidArray, template_parameters="Vector")
    newSolidTensorArray           = PYB11TemplateMethod(newSolidArray, template_parameters="Tensor")
    newSolidSymTensorArray        = PYB11TemplateMethod(newSolidArray, template_parameters="SymTensor")
    newSolidThirdRankTensorArray  = PYB11TemplateMethod(newSolidArray, template_parameters="ThirdRankTensor")
    newSolidFourthRankTensorArray = PYB11TemplateMethod(newSolidArray, template_parameters="FourthRankTensor")
    newSolidFifthRankTensorArray  = PYB11TemplateMethod(newSolidArray, template_parameters="FifthRankTensor")
    newSolidFacetedVolumeArray    = PYB11TemplateMethod(newSolidArray, template_parameters="FacetedVolume")
    newSolidvector_of_intArray    = PYB11TemplateMethod(newSolidArray, template_parameters="std::vector<int>")
    newSolidvector_of_doubleArray = PYB11TemplateMethod(newSolidArray, template_parameters="std::vector<double>")
    newSolidvector_of_VectorArray = PYB11TemplateMethod(newSolidArray, template_parameters="std::vector<Vector>")
    newSolidvector_of_CellFaceFlagArray = PYB11TemplateMethod(newSolidArray, template_parameters="std::vector<CellFaceFlag>")
    newSolidvector_of_CellFaceFlagArray = PYB11TemplateMethod(newSolidArray, template_parameters="std::vector<int>")

    newDEMIntArray              = PYB11TemplateMethod(newDEMArray, template_parameters="int")
    newDEMScalarArray           = PYB11TemplateMethod(newDEMArray, template_parameters="double")
    newDEMVectorArray           = PYB11TemplateMethod(newDEMArray, template_parameters="Vector")
    newDEMTensorArray           = PYB11TemplateMethod(newDEMArray, template_parameters="Tensor")
    newDEMSymTensorArray        = PYB11TemplateMethod(newDEMArray, template_parameters="SymTensor")
    newDEMThirdRankTensorArray  = PYB11TemplateMethod(newDEMArray, template_parameters="ThirdRankTensor")
    newDEMFourthRankTensorArray = PYB11TemplateMethod(newDEMArray, template_parameters="FourthRankTensor")
    newDEMFifthRankTensorArray  = PYB11TemplateMethod(newDEMArray, template_parameters="FifthRankTensor")
    newDEMFacetedVolumeArray    = PYB11TemplateMethod(newDEMArray, template_parameters="FacetedVolume")
    newDEMvector_of_intArray    = PYB11TemplateMethod(newDEMArray, template_parameters="std::vector<int>")
    newDEMvector_of_doubleArray = PYB11TemplateMethod(newDEMArray, template_parameters="std::vector<double>")
    newDEMvector_of_VectorArray = PYB11TemplateMethod(newDEMArray, template_parameters="std::vector<Vector>")
    newDEMvector_of_CellFaceFlagArray = PYB11TemplateMethod(newDEMArray, template_parameters="std::vector<CellFaceFlag>")
    newDEMvector_of_CellFaceFlagArray = PYB11TemplateMethod(newDEMArray, template_parameters="std::vector<int>")

    resizeGlobalIntArray              = PYB11TemplateMethod(resizeGlobalArray, template_parameters="int")
    resizeGlobalScalarArray           = PYB11TemplateMethod(resizeGlobalArray, template_parameters="double")
    resizeGlobalVectorArray           = PYB11TemplateMethod(resizeGlobalArray, template_parameters="Vector")
    resizeGlobalTensorArray           = PYB11TemplateMethod(resizeGlobalArray, template_parameters="Tensor")
    resizeGlobalSymTensorArray        = PYB11TemplateMethod(resizeGlobalArray, template_parameters="SymTensor")
    resizeGlobalThirdRankTensorArray  = PYB11TemplateMethod(resizeGlobalArray, template_parameters="ThirdRankTensor")
    resizeGlobalFourthRankTensorArray = PYB11TemplateMethod(resizeGlobalArray, template_parameters="FourthRankTensor")
    resizeGlobalFifthRankTensorArray  = PYB11TemplateMethod(resizeGlobalArray, template_parameters="FifthRankTensor")
    resizeGlobalFacetedVolumeArray    = PYB11TemplateMethod(resizeGlobalArray, template_parameters="FacetedVolume")
    resizeGlobalvector_of_intArray    = PYB11TemplateMethod(resizeGlobalArray, template_parameters="std::vector<int>")
    resizeGlobalvector_of_doubleArray = PYB11TemplateMethod(resizeGlobalArray, template_parameters="std::vector<double>")
    resizeGlobalvector_of_VectorArray = PYB11TemplateMethod(resizeGlobalArray, template_parameters="std::vector<Vector>")
    resizeGlobalDomainNodeArray       = PYB11TemplateMethod(resizeGlobalArray, template_parameters="DomainNode<%(Dimension)s>")

    resizeFluidIntArray              = PYB11TemplateMethod(resizeFluidArray, template_parameters="int")
    resizeFluidScalarArray           = PYB11TemplateMethod(resizeFluidArray, template_parameters="double")
    resizeFluidVectorArray           = PYB11TemplateMethod(resizeFluidArray, template_parameters="Vector")
    resizeFluidTensorArray           = PYB11TemplateMethod(resizeFluidArray, template_parameters="Tensor")
    resizeFluidSymTensorArray        = PYB11TemplateMethod(resizeFluidArray, template_parameters="SymTensor")
    resizeFluidThirdRankTensorArray  = PYB11TemplateMethod(resizeFluidArray, template_parameters="ThirdRankTensor")
    resizeFluidFourthRankTensorArray = PYB11TemplateMethod(resizeFluidArray, template_parameters="FourthRankTensor")
    resizeFluidFifthRankTensorArray  = PYB11TemplateMethod(resizeFluidArray, template_parameters="FifthRankTensor")
    resizeFluidFacetedVolumeArray    = PYB11TemplateMethod(resizeFluidArray, template_parameters="FacetedVolume")
    resizeFluidvector_of_intArray    = PYB11TemplateMethod(resizeFluidArray, template_parameters="std::vector<int>")
    resizeFluidvector_of_doubleArray = PYB11TemplateMethod(resizeFluidArray, template_parameters="std::vector<double>")
    resizeFluidvector_of_VectorArray = PYB11TemplateMethod(resizeFluidArray, template_parameters="std::vector<Vector>")
    resizeFluidDomainNodeArray       = PYB11TemplateMethod(resizeFluidArray, template_parameters="DomainNode<%(Dimension)s>")

    resizeSolidIntArray              = PYB11TemplateMethod(resizeSolidArray, template_parameters="int")
    resizeSolidScalarArray           = PYB11TemplateMethod(resizeSolidArray, template_parameters="double")
    resizeSolidVectorArray           = PYB11TemplateMethod(resizeSolidArray, template_parameters="Vector")
    resizeSolidTensorArray           = PYB11TemplateMethod(resizeSolidArray, template_parameters="Tensor")
    resizeSolidSymTensorArray        = PYB11TemplateMethod(resizeSolidArray, template_parameters="SymTensor")
    resizeSolidThirdRankTensorArray  = PYB11TemplateMethod(resizeSolidArray, template_parameters="ThirdRankTensor")
    resizeSolidFourthRankTensorArray = PYB11TemplateMethod(resizeSolidArray, template_parameters="FourthRankTensor")
    resizeSolidFifthRankTensorArray  = PYB11TemplateMethod(resizeSolidArray, template_parameters="FifthRankTensor")
    resizeSolidFacetedVolumeArray    = PYB11TemplateMethod(resizeSolidArray, template_parameters="FacetedVolume")
    resizeSolidvector_of_intArray    = PYB11TemplateMethod(resizeSolidArray, template_parameters="std::vector<int>")
    resizeSolidvector_of_doubleArray = PYB11TemplateMethod(resizeSolidArray, template_parameters="std::vector<double>")
    resizeSolidvector_of_VectorArray = PYB11TemplateMethod(resizeSolidArray, template_parameters="std::vector<Vector>")

    resizeDEMIntArray              = PYB11TemplateMethod(resizeDEMArray, template_parameters="int")
    resizeDEMScalarArray           = PYB11TemplateMethod(resizeDEMArray, template_parameters="double")
    resizeDEMVectorArray           = PYB11TemplateMethod(resizeDEMArray, template_parameters="Vector")
    resizeDEMTensorArray           = PYB11TemplateMethod(resizeDEMArray, template_parameters="Tensor")
    resizeDEMSymTensorArray        = PYB11TemplateMethod(resizeDEMArray, template_parameters="SymTensor")
    resizeDEMThirdRankTensorArray  = PYB11TemplateMethod(resizeDEMArray, template_parameters="ThirdRankTensor")
    resizeDEMFourthRankTensorArray = PYB11TemplateMethod(resizeDEMArray, template_parameters="FourthRankTensor")
    resizeDEMFifthRankTensorArray  = PYB11TemplateMethod(resizeDEMArray, template_parameters="FifthRankTensor")
    resizeDEMFacetedVolumeArray    = PYB11TemplateMethod(resizeDEMArray, template_parameters="FacetedVolume")
    resizeDEMvector_of_intArray    = PYB11TemplateMethod(resizeDEMArray, template_parameters="std::vector<int>")
    resizeDEMvector_of_doubleArray = PYB11TemplateMethod(resizeDEMArray, template_parameters="std::vector<double>")
    resizeDEMvector_of_VectorArray = PYB11TemplateMethod(resizeDEMArray, template_parameters="std::vector<Vector>")

    #...........................................................................
    # @PYB11cppname("nodeListPtrs")
    # @PYB11returnpolicy("reference_internal")
    # @PYB11const
    # def nodeLists(self):
    #     return "const std::vector<NodeList<%(Dimension)s>*>&"

    # @PYB11cppname("fluidNodeListPtrs")
    # @PYB11returnpolicy("reference_internal")
    # @PYB11const
    # def fluidNodeLists(self):
    #     return "const std::vector<FluidNodeList<%(Dimension)s>*>&"

    # @PYB11cppname("solidNodeListPtrs")
    # @PYB11returnpolicy("reference_internal")
    # @PYB11const
    # def solidNodeLists(self):
    #     return "const std::vector<SolidNodeList<%(Dimension)s>*>&"

    # @PYB11cppname("DEMNodeListPtrs")
    # @PYB11returnpolicy("reference_internal")
    # @PYB11const
    # def DEMNodeLists(self):
    #     return "const std::vector<DEMNodeList<%(Dimension)s>*>&"

    def setDEMHfieldFromParticleRadius(self, startUniqueIndex = "const int"):
        return "void"

    #...........................................................................
    # Attributes
    nDim = PYB11readonly(static=True, returnpolicy="copy")

    #...........................................................................
    # Properties
    numNodeLists = PYB11property("size_t", "numNodeLists", doc="Number of NodeLists in DataBase")
    numFluidNodeLists = PYB11property("size_t", "numFluidNodeLists", doc="Number of FluidNodeLists in DataBase")
    numSolidNodeLists = PYB11property("size_t", "numSolidNodeLists", doc="Number of SolidNodeLists in DataBase")
    numDEMNodeLists = PYB11property("size_t", "numDEMNodeLists", doc="Number of DEMNodeLists in DataBase")

    numNodes = PYB11property("size_t", "numNodes", doc="Number of nodes in all NodeLists in DataBase")
    numInternalNodes = PYB11property("size_t", "numInternalNodes", doc="Number of internal nodes in all NodeLists in DataBase")
    numGhostNodes = PYB11property("size_t", "numGhostNodes", doc="Number of ghost nodes in all NodeLists in DataBase")

    globalNumNodes = PYB11property("size_t", "globalNumNodes", doc="Number of nodes in all NodeLists in DataBase across all processors")
    globalNumInternalNodes = PYB11property("size_t", "globalNumInternalNodes", doc="Number of internal nodes in all NodeLists in DataBase across all processors")
    globalNumGhostNodes = PYB11property("size_t", "globalNumGhostNodes", doc="Number of ghost nodes in all NodeLists in DataBase across all processors")

    numFluidNodes = PYB11property("size_t", "numFluidNodes", doc="Number of nodes in all NodeLists in DataBase")
    numFluidInternalNodes = PYB11property("size_t", "numFluidInternalNodes", doc="Number of internal nodes in all NodeLists in DataBase")
    numFluidGhostNodes = PYB11property("size_t", "numFluidGhostNodes", doc="Number of ghost nodes in all NodeLists in DataBase")

    globalNumFluidNodes = PYB11property("size_t", "globalNumFluidNodes", doc="Number of nodes in all NodeLists in DataBase across all processors")
    globalNumFluidInternalNodes = PYB11property("size_t", "globalNumFluidInternalNodes", doc="Number of internal nodes in all NodeLists in DataBase across all processors")
    globalNumFluidGhostNodes = PYB11property("size_t", "globalNumFluidGhostNodes", doc="Number of ghost nodes in all NodeLists in DataBase across all processors")
    
    nodeLists = PYB11property("const std::vector<NodeList<%(Dimension)s>*>&", "nodeListPtrs", doc="The set of NodeLists in the DataBase")
    fluidNodeLists = PYB11property("const std::vector<FluidNodeList<%(Dimension)s>*>&", "fluidNodeListPtrs", doc="The set of FluidNodeLists in the DataBase")
    solidNodeLists = PYB11property("const std::vector<SolidNodeList<%(Dimension)s>*>&", "solidNodeListPtrs", doc="The set of SolidNodeLists in the DataBase")
    DEMNodeLists = PYB11property("const std::vector<DEMNodeList<%(Dimension)s>*>&", "DEMNodeListPtrs", doc="The set of NodeLists in the DataBase")
    
    maxKernelExtent = PYB11property("Scalar")
    maxNeighborSearchBuffer = PYB11property("Scalar")

    globalMass = PYB11property("FieldList<%(Dimension)s, Scalar>")
    globalPosition = PYB11property("FieldList<%(Dimension)s, Vector>")
    globalVelocity = PYB11property("FieldList<%(Dimension)s, Vector>")
    globalHfield = PYB11property("FieldList<%(Dimension)s, SymTensor>")
    globalWork = PYB11property("FieldList<%(Dimension)s, Scalar>")

    fluidMass = PYB11property("FieldList<%(Dimension)s, Scalar>")
    fluidPosition = PYB11property("FieldList<%(Dimension)s, Vector>")
    fluidVelocity = PYB11property("FieldList<%(Dimension)s, Vector>")
    fluidMassDensity = PYB11property("FieldList<%(Dimension)s, Scalar>")
    fluidSpecificThermalEnergy = PYB11property("FieldList<%(Dimension)s, Scalar>")
    fluidHfield = PYB11property("FieldList<%(Dimension)s, SymTensor>")
    fluidWork = PYB11property("FieldList<%(Dimension)s, Scalar>")

    solidMass = PYB11property("FieldList<%(Dimension)s, Scalar>")
    solidPosition = PYB11property("FieldList<%(Dimension)s, Vector>")
    solidVelocity = PYB11property("FieldList<%(Dimension)s, Vector>")
    solidMassDensity = PYB11property("FieldList<%(Dimension)s, Scalar>")
    solidSpecificThermalEnergy = PYB11property("FieldList<%(Dimension)s, Scalar>")
    solidHfield = PYB11property("FieldList<%(Dimension)s, SymTensor>")
    solidWork = PYB11property("FieldList<%(Dimension)s, Scalar>")
    solidDeviatoricStress = PYB11property("FieldList<%(Dimension)s, SymTensor>")
    solidPlasticStrain = PYB11property("FieldList<%(Dimension)s, Scalar>")
    solidPlasticStrainRate = PYB11property("FieldList<%(Dimension)s, Scalar>")
    solidDamage = PYB11property("FieldList<%(Dimension)s, SymTensor>")
    solidFragmentIDs = PYB11property("FieldList<%(Dimension)s, int>")

    DEMMass = PYB11property("FieldList<%(Dimension)s, Scalar>")
    DEMPosition = PYB11property("FieldList<%(Dimension)s, Vector>")
    DEMVelocity = PYB11property("FieldList<%(Dimension)s, Vector>")
    DEMParticleRadius = PYB11property("FieldList<%(Dimension)s, Scalar>")
    DEMHfield = PYB11property("FieldList<%(Dimension)s, SymTensor>")
    DEMCompositeParticleIndex = PYB11property("FieldList<%(Dimension)s, int>")
    DEMUniqueIndex = PYB11property("FieldList<%(Dimension)s, size_t>") 

    globalNodeExtent = PYB11property("FieldList<%(Dimension)s, Vector>")
    fluidNodeExtent = PYB11property("FieldList<%(Dimension)s, Vector>")
    solidNodeExtent = PYB11property("FieldList<%(Dimension)s, Vector>")
