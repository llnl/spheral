#-------------------------------------------------------------------------------
# ProbabilisticDamageModel
#-------------------------------------------------------------------------------
from PYB11Generator import *
from DamageModel import *
from RestartMethods import *

@PYB11template("Dimension")
class ProbabilisticDamageModel(DamageModel):
    """ProbabilisticDamageModel

A damage model based on Weibull statistics that uses volume based
probabilities per node to decide when damage starts to accrue.  Should
generate similar results to the classic Benz-Asphaug (Grady-Kipp) model
without generating explicit flaws.  Also appropriate for use with varying
resolution materials."""

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
"""

    def pyinit(self,
               nodeList = "SolidNodeList<%(Dimension)s>&",
               kernel = "const TableKernel<%(Dimension)s>&",
               kWeibull = "const double",
               mWeibull = "const double",
               seed = "const size_t",
               minFlawsPerNode = "const size_t",
               crackGrowthMultiplier = "const double",
               volumeMultiplier = "const double",
               damageCouplingAlgorithm  = "const DamageCouplingAlgorithm",
               strainAlgorithm = "const TensorStrainAlgorithm",
               damageInCompression = "const bool",
               criticalDamageThreshold = "const double",
               mask = "const Field<%(Dimension)s, int>&"):
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
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Compute the derivatives."
        return "void"

    @PYB11virtual
    @PYB11const
    def dt(self,
           dataBase = "const DataBase<%(Dimension)s>&",
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

    @PYB11virtual
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register our state."
        return "void"

    @PYB11virtual
    def registerDerivatives(self,
                            dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Register the derivatives/change fields for updating state."
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
    strainAlgorithm = PYB11property("TensorStrainAlgorithm", "strainAlgorithm",
                                    doc="Choose how the tensor strain is defined to evolve damage")
    damageInCompression = PYB11property("bool", "damageInCompression",
                                        doc="Flag to determine if damage in compression is allowed")
    kWeibull = PYB11property("double", "kWeibull",
                             doc="The coefficient of the Weibull power-law for expected flaw strains")
    mWeibull = PYB11property("double", "mWeibull",
                             doc="The exponent of the Weibull power-law for expected flaw strains")
    volumeMultiplier = PYB11property("double", "volumeMultiplier",
                                     doc="A constant multiplier for the initial node volumes -- meant to mock up higher dimension results in 1 or 2D models.")
    Vmin = PYB11property("double", "Vmin",
                         doc="The minimum initial node volume found")
    Vmax = PYB11property("double", "Vmin",
                         doc="The maximum initial node volume found")
    seed = PYB11property("size_t", "seed",
                         doc="The seed value for the per-node random number generators")
    minFlawsPerNode = PYB11property("size_t", "minFlawsPerNode",
                                    doc="The minimum number of flaws any node will have seeded")
    numFlaws = PYB11property("const Field<%(Dimension)s, int>&", "numFlaws",
                             returnpolicy="reference_internal",
                             doc="The total number of flaws that will be established on each point")
    minFlaw = PYB11property("const Field<%(Dimension)s, Scalar>&", "minFlaw",
                            returnpolicy="reference_internal",
                            doc="The minimum flaw activation strain on each point")
    maxFlaw = PYB11property("const Field<%(Dimension)s, Scalar>&", "maxFlaw",
                            returnpolicy="reference_internal",
                            doc="The maximum flaw activation strain on each point")
    initialVolume = PYB11property("const Field<%(Dimension)s, Scalar>&", "initialVolume",
                                  returnpolicy="reference_internal",
                                  doc="The starting volume for each point")
    youngsModulus = PYB11property("const Field<%(Dimension)s, Scalar>&",
                                  returnpolicy="reference_internal")
    longitudinalSoundSpeed = PYB11property("const Field<%(Dimension)s, Scalar>&",
                                           returnpolicy="reference_internal")
    strain = PYB11property("const Field<%(Dimension)s, SymTensor>&", returnpolicy="reference_internal")
    effectiveStrain = PYB11property("const Field<%(Dimension)s, SymTensor>&", returnpolicy="reference_internal")
    DdamageDt = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    mask = PYB11property("const Field<%(Dimension)s, int>&", "mask", "mask", returnpolicy="reference_internal")
    criticalDamageThreshold = PYB11property("double", "criticalDamageThreshold", "criticalDamageThreshold",
                                            doc="Optional damage threshold, beyond which points do not vote on a timestep")

#-------------------------------------------------------------------------------
# Add the restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, ProbabilisticDamageModel)
