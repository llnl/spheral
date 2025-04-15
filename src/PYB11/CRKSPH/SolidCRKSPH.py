#-------------------------------------------------------------------------------
# SolidCRKSPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from RestartMethods import *
from CRKSPH import *

@PYB11template("Dimension")
@PYB11module("SpheralCRKSPH")
@PYB11dynamic_attr
class SolidCRKSPH(CRKSPH):
    "SolidCRKSPH -- The CRKSPH/ACRKSPH solid material hydrodynamic package for Spheral++."

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ThirdRankTensor = typename %(Dimension)s::ThirdRankTensor;
    using FourthRankTensor = typename %(Dimension)s::FourthRankTensor;
    using FifthRankTensor = typename %(Dimension)s::FifthRankTensor;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
    using PairAccelerationsType = typename SolidCRKSPH<%(Dimension)s>::PairAccelerationsType;
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
               nTensile = "const double",
               damageRelieveRubble = "const bool"):
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
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        """Evaluate the derivatives for the principle hydro variables:
mass density, velocity, and specific thermal energy."""
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

    #...........................................................................
    # Properties
    damageRelieveRubble = PYB11property("bool", "damageRelieveRubble", "damageRelieveRubble",
                                        doc="Control whether allow damaged material to have stress relieved.")
    DdeviatoricStressDt = PYB11property("const FieldList<%(Dimension)s, SymTensor>&", "DdeviatoricStressDt", returnpolicy="reference_internal")
    bulkModulus = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "bulkModulus", returnpolicy="reference_internal")
    shearModulus = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "shearModulus", returnpolicy="reference_internal")
    yieldStrength = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "yieldStrength", returnpolicy="reference_internal")
    plasticStrain0 = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "plasticStrain0", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, SolidCRKSPH)

