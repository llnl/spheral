#-------------------------------------------------------------------------------
# CRKSPHBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from GenericHydro import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCRKSPH")
@PYB11dynamic_attr
class CRKSPHBase(GenericHydro):
    "CRKSPHBase -- Base class for the CRKSPH/ACRKSPH hydrodynamic packages"

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
"""

    def pyinit(self,
               dataBase = "DataBase<%(Dimension)s>&",
               Q = "ArtificialViscosityHandle<%(Dimension)s>&",
               order = "const RKOrder",
               cfl = "const double",
               useVelocityMagnitudeForDt = "const bool",
               compatibleEnergyEvolution = "const bool",
               evolveTotalEnergy = "const bool",
               XSPH = "const bool",
               densityUpdate = "const MassDensityType",
               epsTensile = "const double",
               nTensile = "const double"):
        "Constructor"

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
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state Hydro expects to use and evolve."
        return "void"

    @PYB11virtual
    def registerDerivatives(self,
                            dataBase = "DataBase<%(Dimension)s>&",
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
    def finalizeDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Finalize the derivatives."
        return "void"

    @PYB11virtual
    def applyGhostBoundaries(self,
                             state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    def enforceBoundaries(self,
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Enforce boundary conditions for the physics specific fields."
        return "void"

    @PYB11virtual
    @PYB11const
    def requireReproducingKernels(self):
        "CRK overrides and provides this method"
        return "std::set<RKOrder>"

    #...........................................................................
    # Properties
    correctionOrder = PYB11property("RKOrder", "correctionOrder", "correctionOrder",
                                    doc="Flag to choose CRK Correction Order")
    densityUpdate = PYB11property("MassDensityType", "densityUpdate", "densityUpdate",
                                  doc="Flag to choose whether we want to sum for density, or integrate the continuity equation.")
    compatibleEnergyEvolution = PYB11property("bool", "compatibleEnergyEvolution", "compatibleEnergyEvolution",
                                              doc="Flag to determine if we're using the total energy conserving compatible energy evolution scheme.")
    evolveTotalEnergy = PYB11property("bool", "evolveTotalEnergy", "evolveTotalEnergy",
                                      doc="Flag controlling if we evolve total or specific energy.")
    XSPH = PYB11property("bool", "XSPH", "XSPH",
                         doc="Flag to determine if we're using the XSPH algorithm.")
    epsilonTensile = PYB11property("Scalar", "epsilonTensile", "epsilonTensile",
                                   doc="Parameters for the tensile correction force at small scales.")
    nTensile = PYB11property("Scalar", "nTensile", "nTensile",
                                   doc="Parameters for the tensile correction force at small scales.")

    timeStepMask = PYB11property("const FieldList<%(Dimension)s, int>&", "timeStepMask", returnpolicy="reference_internal")
    pressure = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "pressure", returnpolicy="reference_internal")
    soundSpeed = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "soundSpeed", returnpolicy="reference_internal")
    entropy = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "entropy", returnpolicy="reference_internal")
    XSPHDeltaV = PYB11property("const FieldList<%(Dimension)s, Vector>&", "XSPHDeltaV", returnpolicy="reference_internal")

    DxDt = PYB11property("const FieldList<%(Dimension)s, Vector>&", "DxDt", returnpolicy="reference_internal")
    DvDt = PYB11property("const FieldList<%(Dimension)s, Vector>&", "DvDt", returnpolicy="reference_internal")
    DmassDensityDt = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "DmassDensityDt", returnpolicy="reference_internal")
    DspecificThermalEnergyDt = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "DspecificThermalEnergyDt", returnpolicy="reference_internal")
    DvDx = PYB11property("const FieldList<%(Dimension)s, Tensor>&", "DvDx", returnpolicy="reference_internal")
    internalDvDx = PYB11property("const FieldList<%(Dimension)s, Tensor>&", "internalDvDx", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, CRKSPHBase)

