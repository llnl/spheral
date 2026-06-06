#-------------------------------------------------------------------------------
# GenericRiemannHydro
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralGSPH")
@PYB11dynamic_attr
class GenericRiemannHydro(Physics):

    PYB11typedefs = """
  using Scalar = typename %(Dimension)s::Scalar;
  using Vector = typename %(Dimension)s::Vector;
  using Tensor = typename %(Dimension)s::Tensor;
  using SymTensor = typename %(Dimension)s::SymTensor;
  using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
  using PairAccelerationsType = typename GenericRiemannHydro<%(Dimension)s>::PairAccelerationsType;
  using PairWorkType = typename GenericRiemannHydro<%(Dimension)s>::PairWorkType;
  using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
"""
    
    def pyinit(dataBase = "DataBase<%(Dimension)s>&",
               riemannSolver = "RiemannSolverBase<%(Dimension)s>&",
               W = "const TableKernel<%(Dimension)s>&",
               epsDiffusionCoeff = "const Scalar",
               cfl = "const double",
               useVelocityMagnitudeForDt = "const bool",
               compatibleEnergyEvolution = "const bool",
               evolveTotalEnergy = "const bool",
               XSPH = "const bool",
               correctVelocityGradient = "const bool",
               gradType = "const GradientType",
               densityUpdate = "const MassDensityType",
               epsTensile = "const double",
               nTensile = "const double",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "GenericRiemannHydro constructor"

    #...........................................................................
    # Virtual methods

    @PYB11virtual
    @PYB11const
    def dt(dataBase = "const DataBase<%(Dimension)s>&", 
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

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
        "Register the state Hydro expects to use and evolve."
        return "void"


    @PYB11virtual
    def preStepInitialize(self,
                          dataBase = "const DataBase<%(Dimension)s>&", 
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Optional hook to be called at the beginning of a time step."
        return "void"

    @PYB11virtual
    def initialize(time = "const Scalar",
                   dt = "const Scalar",
                   dataBase = "const DataBase<%(Dimension)s>&",
                   state = "State<%(Dimension)s>&",
                   derivs = "StateDerivatives<%(Dimension)s>&"):
        "Initialize the Hydro before we start a derivative evaluation."
        return "bool"
                       
#     @PYB11virtual
#     @PYB11const
#     def evaluateDerivatives(time = "const Scalar",
#                             dt = "const Scalar",
#                             dataBase = "const DataBase<%(Dimension)s>&",
#                             state = "const State<%(Dimension)s>&",
#                             derivs = "StateDerivatives<%(Dimension)s>&"):
#         """Evaluate the derivatives for the principle hydro 
# mass density, velocity, and specific thermal energy."""
#         return "void"

    @PYB11virtual
    @PYB11const
    def finalizeDerivatives(time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Finalize the derivatives."
        return "void"

    @PYB11virtual
    def applyGhostBoundaries(state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    def enforceBoundaries(state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Enforce boundary conditions for the physics specific fields."
        return "void"


    # #...........................................................................
    # # Protected methods from GenericHydro
    # @PYB11protected
    # @PYB11const
    # def updateMasterNeighborStats(self, numMaster="int"):
    #     return "void"

    # @PYB11protected
    # @PYB11const
    # def updateCoarseNeighborStats(self, numCoarse="int"):
    #     return "void"

    # @PYB11protected
    # @PYB11const
    # def updateRefineNeighborStats(self, numRefine="int"):
    #     return "void"

    # @PYB11protected
    # @PYB11const
    # def updateActualNeighborStats(self, numActual="int"):
    #     return "void"


    #...........................................................................
    # Properties
    cfl = PYB11property("Scalar", "cfl", "cfl", doc="The Courant-Friedrichs-Lewy timestep limit multiplier")
    specificThermalEnergyDiffusionCoefficient = PYB11property("Scalar", "specificThermalEnergyDiffusionCoefficient", "specificThermalEnergyDiffusionCoefficient", 
                                          doc="coefficient used to diffuse specificThermalEnergy amongst like nodes.")
    riemannSolver = PYB11property("RiemannSolverBase<%(Dimension)s>&", "riemannSolver", doc="The object defining the interface state construction.")
    kernel = PYB11property("const TableKernel<%(Dimension)s>&", "kernel", doc="The interpolation kernel")
    gradientType = PYB11property("GradientType", "gradientType", "gradientType",
                                 doc="Enum to selecting different gradients we can use")
    densityUpdate = PYB11property("MassDensityType", "densityUpdate", "densityUpdate",
                                  doc="Flag to choose whether we want to sum for density, or integrate the continuity equation.")
    compatibleEnergyEvolution = PYB11property("bool", "compatibleEnergyEvolution", "compatibleEnergyEvolution",
                                              doc="Flag to determine if we're using the total energy conserving compatible energy evolution scheme.")
    evolveTotalEnergy = PYB11property("bool", "evolveTotalEnergy", "evolveTotalEnergy",
                                      doc="Flag controlling if we evolve total or specific energy.")
    XSPH = PYB11property("bool", "XSPH", "XSPH",
                         doc="Flag to determine if we're using the XSPH algorithm.")
    correctVelocityGradient = PYB11property("bool", "correctVelocityGradient", "correctVelocityGradient",
                                            doc="Flag to determine if we're applying the linear correction for the velocity gradient.")
    epsilonTensile = PYB11property("double", "epsilonTensile", "epsilonTensile",
                                   doc="Parameters for the tensile correction force at small scales.")
    nTensile = PYB11property("double", "nTensile", "nTensile",
                             doc="Parameters for the tensile correction force at small scales.")
    xmin = PYB11property("const Vector&", "xmin", "xmin",
                         returnpolicy="reference_internal",
                         doc="Optional minimum coordinate for bounding box for use generating the mesh for the Voronoi mass density update.")
    xmax = PYB11property("const Vector&", "xmax", "xmax",
                         returnpolicy="reference_internal",
                         doc="Optional maximum coordinate for bounding box for use generating the mesh for the Voronoi mass density update.")

    timeStepMask =                 PYB11property("const FieldList<%(Dimension)s, int>&",      "timeStepMask",         returnpolicy="reference_internal")
    pressure =                     PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "pressure",             returnpolicy="reference_internal")
    volume =                       PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "volume",             returnpolicy="reference_internal")
    soundSpeed =                   PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "soundSpeed",           returnpolicy="reference_internal")
    normalization =                PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "normalization",        returnpolicy="reference_internal")
    XSPHWeightSum =                PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "XSPHWeightSum",        returnpolicy="reference_internal")
    XSPHDeltaV =                   PYB11property("const FieldList<%(Dimension)s, Vector>&",   "XSPHDeltaV",           returnpolicy="reference_internal")
    M =                            PYB11property("const FieldList<%(Dimension)s, Tensor>&",   "M",                    returnpolicy="reference_internal")
    DxDt =                         PYB11property("const FieldList<%(Dimension)s, Vector>&",   "DxDt",                 returnpolicy="reference_internal")
    DvDt =                         PYB11property("const FieldList<%(Dimension)s, Vector>&",   "DvDt",                 returnpolicy="reference_internal")
    DspecificThermalEnergyDt =     PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "DspecificThermalEnergyDt", returnpolicy="reference_internal")
    DvDx =                         PYB11property("const FieldList<%(Dimension)s, Tensor>&",   "DvDx",                 returnpolicy="reference_internal")
    
    pairAccelerations = PYB11property("const PairAccelerationsType&", "pairAccelerations", returnpolicy="reference_internal")
    pairDepsDt = PYB11property("const PairWorkType&", "pairDepsDt", returnpolicy="reference_internal")
    riemannDpDx = PYB11property("const FieldList<%(Dimension)s, Vector>&","riemannDpDx",returnpolicy="reference_internal")
    newRiemannDpDx = PYB11property("const FieldList<%(Dimension)s, Vector>&","newRiemannDpDx",returnpolicy="reference_internal")
    riemannDvDx = PYB11property("const FieldList<%(Dimension)s, Tensor>&","riemannDvDx",returnpolicy="reference_internal")
    newRiemannDvDx = PYB11property("const FieldList<%(Dimension)s, Tensor>&","newRiemannDvDx",returnpolicy="reference_internal")
    
    # useVelocityMagnitudeForDt = PYB11property("bool", "useVelocityMagnitudeForDt", "useVelocityMagnitudeForDt", doc="Should the pointwise velocity magnitude be used to limit the timestep?")
    # minMasterNeighbor = PYB11property("int", "minMasterNeighbor", doc="minimum number of master neighbors found")
    # maxMasterNeighbor = PYB11property("int", "maxMasterNeighbor", doc="maximum number of master neighbors found")
    # averageMasterNeighbor = PYB11property("double", "averageMasterNeighbor", doc="average number of master neighbors found")
    # minCoarseNeighbor = PYB11property("int", "minCoarseNeighbor", doc="minimum number of coarse neighbors found")
    # maxCoarseNeighbor = PYB11property("int", "maxCoarseNeighbor", doc="maximum number of coarse neighbors found")
    # averageCoarseNeighbor = PYB11property("double", "averageCoarseNeighbor", doc="average number of coarse neighbors found")
    # minRefineNeighbor = PYB11property("int", "minRefineNeighbor", doc="minimum number of refine neighbors found")
    # maxRefineNeighbor = PYB11property("int", "maxRefineNeighbor", doc="maximum number of refine neighbors found")
    # averageRefineNeighbor = PYB11property("double", "averageRefineNeighbor", doc="average number of refine neighbors found")
    # minActualNeighbor = PYB11property("int", "minActualNeighbor", doc="minimum number of actual neighbors found")
    # maxActualNeighbor = PYB11property("int", "maxActualNeighbor", doc="maximum number of actual neighbors found")
    # averageActualNeighbor = PYB11property("double", "averageActualNeighbor", doc="average number of actual neighbors found")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, GenericRiemannHydro)
