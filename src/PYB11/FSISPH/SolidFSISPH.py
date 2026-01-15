#-------------------------------------------------------------------------------
# SolidFSISPH
#-------------------------------------------------------------------------------
from PYB11Generator import *
from GenericHydro import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralFSISPH")
@PYB11dynamic_attr
class SolidFSISPH(GenericHydro):
    "SolidFSISPH -- SolidSPHHydro modified for large density discontinuities"

    PYB11typedefs = """
  using Scalar = typename %(Dimension)s::Scalar;
  using Vector = typename %(Dimension)s::Vector;
  using Tensor = typename %(Dimension)s::Tensor;
  using SymTensor = typename %(Dimension)s::SymTensor;
  using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
  using PairAccelerationsType = typename SolidFSISPH<%(Dimension)s>::PairAccelerationsType;
  using PairWorkType = typename SolidFSISPH<%(Dimension)s>::PairWorkType;
  using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
"""
    
    def pyinit(dataBase = "DataBase<%(Dimension)s>&",
               Q = "ArtificialViscosity<%(Dimension)s>&",
               slides = "SlideSurface<%(Dimension)s>&",
               W = "const TableKernel<%(Dimension)s>&",
               cfl = "const double",
               surfaceForceCoefficient = "const double",
               densityStabilizationCoefficient = "const double",
               specificThermalEnergyDiffusionCoefficient = "const double",
               xsphCoefficient = "const double",
               interfaceMethod = "const InterfaceMethod",
               kernelAveragingMethod = "const KernelAveragingMethod",
               sumDensityNodeLists = "std::vector<int>",
               useVelocityMagnitudeForDt = "const bool",
               compatibleEnergyEvolution = "const bool",
               evolveTotalEnergy = "const bool",
               linearCorrectGradients = "const bool",
               decoupleDamagedMaterial = "const bool",
               interfacePmin = "const double",
               interfaceNeighborAngleThreshold = "const double ",
               densityUpdate = "const FSIMassDensityMethod",
               epsTensile = "const double",
               nTensile = "const double",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "SolidFSISPH constructor"

    #...........................................................................
    # Virtual methods
    
    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        """Evaluate the derivatives for the principle hydro 
mass density, velocity, and specific thermal energy."""
        return "void"

    @PYB11virtual
    def initializeProblemStartupDependencies(self,
                                             dataBase = "DataBase<%(Dimension)s>&",
                                             state = "State<%(Dimension)s>&",
                                             derivs = "StateDerivatives<%(Dimension)s>&"):
        """A second optional method to be called on startup, after Physics::initializeProblemStartup has
been called.
One use for this hook is to fill in dependendent state using the State object, such as
temperature or pressure."""
        return "void"

    @PYB11virtual
    def registerState(dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "register the surface normals"
        return "void"


    @PYB11virtual
    def registerDerivatives(dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "non-op place filler"
        return "void"

    #...........................................................................
    # Properties
    kernel = PYB11property("const TableKernel<%(Dimension)s>&", "kernel", doc="The interpolation kernel")
    slideSurfaces = PYB11property("SlideSurface<%(Dimension)s>&", "slideSurface", doc="The slide surface object")
    
    densityUpdate = PYB11property("FSIMassDensityMethod", "densityUpdate", "densityUpdate", doc="Flag to choose whether we want to sum for density, or integrate the continuity equation.")
    interfaceMethod = PYB11property("InterfaceMethod", "interfaceMethod", "interfaceMethod",doc="Flag to select how we want construct material interfaces")
    kernelAveragingMethod = PYB11property("KernelAveragingMethod", "kernelAveragingMethod", "kernelAveragingMethod",doc="Flag to select our kernel type")

    planeStrain = PYB11property("bool", "planeStrain", "planeStrain",doc="use plane strain approach for 1D or 2D problems.")
    decoupleDamagedMaterial = PYB11property("bool", "decoupleDamagedMaterial", "decoupleDamagedMaterial",doc="turn off decoupling algo for same-material damaged interactions.")
    compatibleEnergyEvolution = PYB11property("bool", "compatibleEnergyEvolution", "compatibleEnergyEvolution",doc="Flag to determine if we're using the total energy conserving compatible energy evolution scheme.")
    evolveTotalEnergy = PYB11property("bool", "evolveTotalEnergy", "evolveTotalEnergy",doc="Flag controlling if we evolve total or specific energy.")
    linearCorrectGradients = PYB11property("bool", "linearCorrectGradients", "linearCorrectGradients",doc="Flag to determine if we're applying the linear correction for the velocity gradient.")
    sumDensityNodeLists = PYB11property("std::vector<int>", "sumDensityNodeLists", "sumDensityNodeLists",doc="control if rigorous density sum is applied to individual node lists.")
    
    surfaceForceCoefficient = PYB11property("double", "surfaceForceCoefficient", "surfaceForceCoefficient",doc="additional force between different materials ala Monaghan 2013.")
    densityStabilizationCoefficient = PYB11property("double", "densityStabilizationCoefficient", "densityStabilizationCoefficient",doc="coefficient used to adjust velocity gradient to prevent unstable rho.")
    specificThermalEnergyDiffusionCoefficient = PYB11property("double", "specificThermalEnergyDiffusionCoefficient", "specificThermalEnergyDiffusionCoefficient",doc="coefficient used to diffuse specificThermalEnergy amongst like nodes.")
    xsphCoefficient = PYB11property("double", "xsphCoefficient", "xsphCoefficient",doc="coefficient to dial magnitude of xsph.")
    interfacePmin = PYB11property("double", "interfacePmin", "interfacePmin",doc="minimum pressure allowed across a material face.")
    interfaceNeighborAngleThreshold = PYB11property("double", "interfaceNeighborAngleThreshold", "interfaceNeighborAngleThreshold",doc="parameter controling interface and free surface detection.")

    epsilonTensile = PYB11property("double", "epsilonTensile", "epsilonTensile",doc="Parameters for the tensile correction force at small scales.")
    nTensile = PYB11property("double", "nTensile", "nTensile", doc="Parameters for the tensile correction force at small scales.")
    xmin = PYB11property("const Vector&", "xmin", "xmin",returnpolicy="reference_internal",doc="Optional minimum coordinate for bounding box for use generating the mesh for the Voronoi mass density update.")
    xmax = PYB11property("const Vector&", "xmax", "xmax",returnpolicy="reference_internal",doc="Optional maximum coordinate for bounding box for use generating the mesh for the Voronoi mass density update.")
    
    pairAccelerations = PYB11property("const PairAccelerationsType&", "pairAccelerations", returnpolicy="reference_internal")
    pairDepsDt = PYB11property("const PairWorkType&", "pairDepsDt", returnpolicy="reference_internal")

    timeStepMask =                 PYB11property("const FieldList<%(Dimension)s, int>&",      "timeStepMask",         returnpolicy="reference_internal")
    pressure =                     PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "pressure",             returnpolicy="reference_internal")
    damagedPressure =              PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "damagedPressure",      returnpolicy="reference_internal")
    soundSpeed =                   PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "soundSpeed",           returnpolicy="reference_internal")
    volume =                       PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "volume",               returnpolicy="reference_internal")
    bulkModulus =                  PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "bulkModulus",         returnpolicy="reference_internal")
    shearModulus =                 PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "shearModulus",        returnpolicy="reference_internal")
    yieldStrength =                PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "yieldStrength",       returnpolicy="reference_internal")
    plasticStrain0 =               PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "plasticStrain0",      returnpolicy="reference_internal")
    #inverseEquivalentDeviatoricStress = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "inverseEquivalentDeviatoricStress",returnpolicy="reference_internal")
    DxDt =                         PYB11property("const FieldList<%(Dimension)s, Vector>&",   "DxDt",                 returnpolicy="reference_internal")
    XSPHDeltaV =                   PYB11property("const FieldList<%(Dimension)s, Vector>&",   "XSPHDeltaV",           returnpolicy="reference_internal")
    XSPHWeightSum =                PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "XSPHWeightSum",        returnpolicy="reference_internal")
    DvDt =                         PYB11property("const FieldList<%(Dimension)s, Vector>&",   "DvDt",                 returnpolicy="reference_internal")
    DmassDensityDt =               PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "DmassDensityDt",       returnpolicy="reference_internal")
    DspecificThermalEnergyDt =     PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "DspecificThermalEnergyDt", returnpolicy="reference_internal")
    DdeviatoricStressDt =          PYB11property("const FieldList<%(Dimension)s, SymTensor>&","DdeviatoricStressDt",  returnpolicy="reference_internal")
    DepsDx =                       PYB11property("const FieldList<%(Dimension)s, Vector>&",   "DepsDx",               returnpolicy="reference_internal")
    DPDx =                         PYB11property("const FieldList<%(Dimension)s, Vector>&",   "DPDx",                 returnpolicy="reference_internal")
    DvDx =                         PYB11property("const FieldList<%(Dimension)s, Tensor>&",   "DvDx",                 returnpolicy="reference_internal")
    internalDvDx =                 PYB11property("const FieldList<%(Dimension)s, Tensor>&",   "internalDvDx",         returnpolicy="reference_internal")
    M =                            PYB11property("const FieldList<%(Dimension)s, Tensor>&",   "M",                    returnpolicy="reference_internal")
    localM =                       PYB11property("const FieldList<%(Dimension)s, Tensor>&",   "localM",               returnpolicy="reference_internal")
    normalization =                PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "normalization",        returnpolicy="reference_internal")
    interfaceFraction =                PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "interfaceFraction",    returnpolicy="reference_internal")
    interfaceFlags =                   PYB11property("const FieldList<%(Dimension)s, int>&",      "interfaceFlags",       returnpolicy="reference_internal")
    interfaceAreaVectors =             PYB11property("const FieldList<%(Dimension)s, Vector>&",   "interfaceAreaVectors", returnpolicy="reference_internal")
    interfaceNormals =                 PYB11property("const FieldList<%(Dimension)s, Vector>&",   "interfaceNormals",     returnpolicy="reference_internal")
    interfaceSmoothness =              PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "interfaceSmoothness",  returnpolicy="reference_internal")
    newInterfaceAreaVectors =          PYB11property("const FieldList<%(Dimension)s, Vector>&",   "newInterfaceAreaVectors", returnpolicy="reference_internal")
    newInterfaceNormals =              PYB11property("const FieldList<%(Dimension)s, Vector>&",   "newInterfaceNormals",  returnpolicy="reference_internal")
    interfaceSmoothnessNormalization = PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "interfaceSmoothnessNormalization", returnpolicy="reference_internal")
    newInterfaceSmoothness =           PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "newInterfaceSmoothness",returnpolicy="reference_internal")
    interfaceAngles =                  PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "interfaceAngles",       returnpolicy="reference_internal")
    
#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, SolidFSISPH)
