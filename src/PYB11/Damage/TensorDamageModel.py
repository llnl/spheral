#-------------------------------------------------------------------------------
# TensorDamageModel
#-------------------------------------------------------------------------------
from PYB11Generator import *
from DamageModel import *
from RestartMethods import *

@PYB11template("Dimension")
class TensorDamageModel(DamageModel):
    """TensorDamageModel -- Base class for the tensor damage physics models.
This class does not know how to seed the flaw distribution -- that is 
required of descendant classes."""

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;

    using FlawStorageType = Field<%(Dimension)s, std::vector<double> >;
"""

    def pyinit(self,
               nodeList = "SolidNodeList<%(Dimension)s>&",
               strainAlgorithm = "const TensorStrainAlgorithm",
               damageCouplingAlgorithm  = "const DamageCouplingAlgorithm",
               kernel = "const TableKernel<%(Dimension)s>&",
               crackGrowthMultiplier = "const double",
               criticalDamageThreshold = "const double",
               damageInCompression = "const bool",
               flaws = "const FlawStorageType&"):
        "Constructor"

    def pyinit1(self,
                nodeList = "SolidNodeList<%(Dimension)s>&",
                strainAlgorithm = "const TensorStrainAlgorithm",
                damageCouplingAlgorithm  = "const DamageCouplingAlgorithm",
                kernel = "const TableKernel<%(Dimension)s>&",
                crackGrowthMultiplier = "const double",
                criticalDamageThreshold = "const double",
                damageInCompression = "const bool"):
        "Constructor"
    #...........................................................................
    # Virtual methods
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

    #...........................................................................
    # Methods
    def cullToWeakestFlaws(sefl):
        "Optional method to cull the set of flaws to the single weakest one on each point."
        return "void"

    @PYB11const
    def flawsForNode(self, index="const size_t"):
        "Get the set of flaw activation energies for the given node index."
        return "const std::vector<double>"

    #...........................................................................
    # Properties
    youngsModulus = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    longitudinalSoundSpeed = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    flaws = PYB11property("const FlawStorageType&", "flaws", "flaws", returnpolicy="reference_internal",
                          doc="The raw set of flaw activation strains per point")
    sumActivationEnergiesPerNode = PYB11property("Field<%(Dimension)s, Scalar>", 
                                                 doc="Compute a Field with the sum of the activation energies per node.")
    numFlawsPerNode = PYB11property("Field<%(Dimension)s, Scalar>",
                                    doc="Compute a Field with the number of flaws per node.")
    strain = PYB11property("const Field<%(Dimension)s, SymTensor>&", returnpolicy="reference_internal")
    effectiveStrain = PYB11property("const Field<%(Dimension)s, SymTensor>&", returnpolicy="reference_internal")
    DdamageDt = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")

    strainAlgorithm = PYB11property("TensorStrainAlgorithm")
    damageInCompression = PYB11property("bool", "damageInCompression", "damageInCompression",
                                        doc="Flag to determine if damage in compression is allowed.")
    criticalDamageThreshold = PYB11property("double", "criticalDamageThreshold", "criticalDamageThreshold",
                                            doc="The critical damage threshold for not setting the time step.")

#-------------------------------------------------------------------------------
# Add the restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, TensorDamageModel)
