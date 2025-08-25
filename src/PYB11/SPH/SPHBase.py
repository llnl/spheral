#-------------------------------------------------------------------------------
# SPHBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from GenericHydro import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralSPH")
@PYB11dynamic_attr
class SPHBase(GenericHydro):

    PYB11typedefs = """
  using Scalar = typename %(Dimension)s::Scalar;
  using Vector = typename %(Dimension)s::Vector;
  using Tensor = typename %(Dimension)s::Tensor;
  using SymTensor = typename %(Dimension)s::SymTensor;
  using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
  using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
"""
    
    def pyinit(dataBase = "DataBase<%(Dimension)s>&",
               Q = "ArtificialViscosity<%(Dimension)s>&",
               W = "const TableKernel<%(Dimension)s>&",
               WPi = "const TableKernel<%(Dimension)s>&",
               cfl = "const double",
               useVelocityMagnitudeForDt = "const bool",
               compatibleEnergyEvolution = "const bool",
               evolveTotalEnergy = "const bool",
               gradhCorrection = "const bool",
               XSPH = "const bool",
               correctVelocityGradient = "const bool",
               sumMassDensityOverAllNodeLists = "const bool",
               densityUpdate = "const MassDensityType",
               epsTensile = "const double",
               nTensile = "const double",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "SPHBase constructor"

    #...........................................................................
    # Virtual methods
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
    def registerDerivatives(dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Register the derivatives/change fields for updating state."
        return "void"

    @PYB11virtual
    def preStepInitialize(self,
                          dataBase = "const DataBase<%(Dimension)s>&", 
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Optional hook to be called at the beginning of a time step."
        return "void"

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
    def postStateUpdate(time = "const Scalar",
                        dt = "const Scalar",
                        dataBase = "const DataBase<%(Dimension)s>&",
                        state = "State<%(Dimension)s>&",
                        derivs = "StateDerivatives<%(Dimension)s>&"):
        "Post-state update. For SPH this is where we recompute the omega (grad h) corrections."
        return "bool"
               
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

    #...........................................................................
    # Methods
    @PYB11const
    def updateVolume(state = "State<%(Dimension)s>&",
                     boundaries = "const bool"):
        """A method to fill in the volume in the State, optionally enforcing
boundary conditions."""
        return "void"

    #...........................................................................
    # Properties
    kernel = PYB11property("const TableKernel<%(Dimension)s>&", "kernel", doc="The interpolation kernel")
    PiKernel = PYB11property("const TableKernel<%(Dimension)s>&", "PiKernel", doc="The interpolation kernel for the artificial viscosity")
    densityUpdate = PYB11property("MassDensityType", "densityUpdate", "densityUpdate",
                                  doc="Flag to choose whether we want to sum for density, or integrate the continuity equation.")
    compatibleEnergyEvolution = PYB11property("bool", "compatibleEnergyEvolution", "compatibleEnergyEvolution",
                                              doc="Flag to determine if we're using the total energy conserving compatible energy evolution scheme.")
    evolveTotalEnergy = PYB11property("bool", "evolveTotalEnergy", "evolveTotalEnergy",
                                      doc="Flag controlling if we evolve total or specific energy.")
    gradhCorrection = PYB11property("bool", "gradhCorrection", "gradhCorrection",
                                    doc="Flag to determine if we're using the grad h correction.")
    XSPH = PYB11property("bool", "XSPH", "XSPH",
                         doc="Flag to determine if we're using the XSPH algorithm.")
    correctVelocityGradient = PYB11property("bool", "correctVelocityGradient", "correctVelocityGradient",
                                            doc="Flag to determine if we're applying the linear correction for the velocity gradient.")
    sumMassDensityOverAllNodeLists = PYB11property("bool", "sumMassDensityOverAllNodeLists", "sumMassDensityOverAllNodeLists",
                                                   doc="Flag to determine if the sum density definition extends over neighbor NodeLists.")
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
    soundSpeed =                   PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "soundSpeed",           returnpolicy="reference_internal")
    volume =                       PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "volume",               returnpolicy="reference_internal")
    omegaGradh =                   PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "omegaGradh",           returnpolicy="reference_internal")
    entropy =                      PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "entropy",              returnpolicy="reference_internal")
    massDensityCorrection =        PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "massDensityCorrection",returnpolicy="reference_internal")
    massDensitySum =               PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "massDensitySum",       returnpolicy="reference_internal")
    normalization =                PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "normalization",        returnpolicy="reference_internal")
    XSPHWeightSum =                PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "XSPHWeightSum",        returnpolicy="reference_internal")
    XSPHDeltaV =                   PYB11property("const FieldList<%(Dimension)s, Vector>&",   "XSPHDeltaV",           returnpolicy="reference_internal")
    M =                            PYB11property("const FieldList<%(Dimension)s, Tensor>&",   "M",                    returnpolicy="reference_internal")
    localM =                       PYB11property("const FieldList<%(Dimension)s, Tensor>&",   "localM",               returnpolicy="reference_internal")
    DxDt =                         PYB11property("const FieldList<%(Dimension)s, Vector>&",   "DxDt",                 returnpolicy="reference_internal")
    DvDt =                         PYB11property("const FieldList<%(Dimension)s, Vector>&",   "DvDt",                 returnpolicy="reference_internal")
    DmassDensityDt =               PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "DmassDensityDt",       returnpolicy="reference_internal")
    DspecificThermalEnergyDt =     PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "DspecificThermalEnergyDt", returnpolicy="reference_internal")
    DvDx =                         PYB11property("const FieldList<%(Dimension)s, Tensor>&",   "DvDx",                 returnpolicy="reference_internal")
    internalDvDx =                 PYB11property("const FieldList<%(Dimension)s, Tensor>&",   "internalDvDx",         returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, SPHBase)
