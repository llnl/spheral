#-------------------------------------------------------------------------------
# DamageModel
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from RestartMethods import *

@PYB11template("Dimension")
class DamageModel(Physics):
    """DamageModel -- Base class for the damage physics models.
This class just provides the basic interface for damage models, and does 
not fill out the complete physics package interface."""

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
               W = "const TableKernel<%(Dimension)s>&",
               crackGrowthMultiplier = "const double",
               damageCouplingAlgorithm  = "const DamageCouplingAlgorithm"):
        "Constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def computeScalarDDDt(self,
                          dataBase = "const DataBase<%(Dimension)s>&",
                          state = "const State<%(Dimension)s>&",
                          time = "const Scalar",
                          dt = "const Scalar",
                          DDDt = "Field<%(Dimension)s, Scalar>&"):
        "Compute the generic Grady-Kipp (ala Benz-Asphaug) scalar damage time derivative."
        return "void"

    @PYB11virtual 
    def initialize(self,
                   time = "const Scalar",
                   dt = "const Scalar",
                   dataBase = "const DataBase<%(Dimension)s>&",
                   state = "State<%(Dimension)s>&",
                   derivs = "StateDerivatives<%(Dimension)s>&"):
        "Initialize before computing the derivatives."
        return "bool"

    @PYB11virtual
    def finalize(self,
                 time = "const Scalar", 
                 dt = "const Scalar",
                 dataBase = "DataBase<%(Dimension)s>&", 
                 state = "State<%(Dimension)s>&",
                 derivs = "StateDerivatives<%(Dimension)s>&"):
        "Similarly packages might want a hook to do some post-step finalizations.  Really we should rename this post-step finalize."
        return "void"

    @PYB11virtual
    @PYB11const
    def requireGhostConnectivity(self):
        "Some physics algorithms require ghost connectivity to be constructed."
        return "bool"

    @PYB11virtual
    @PYB11const
    def requireIntersectionConnectivity(self):
        "Some physics algorithms require intersection connectivity to be constructed."
        return "bool"

    #...........................................................................
    # Properties
    nodeList = PYB11property("const SolidNodeList<%(Dimension)s>&", returnpolicy="reference_internal",
                             doc="Access the SolidNodeList we're damaging.")
    kernel = PYB11property("const TableKernel<%(Dimension)s>&", returnpolicy="reference_internal",
                           doc="Access the kernel.")
    crackGrowthMultiplier = PYB11property("double")
    damageCouplingAlgorithm = PYB11property("DamageCouplingAlgorithm", "damageCouplingAlgorithm",
                                            doc="The choice for coupling with damaged nodes.")
    nodeCoupling = PYB11property("const NodeCoupling&", "nodeCoupling",
                                 doc="The NodeCoupling implementation")
    excludeNodes = PYB11property("std::vector<int>", "excludeNodes", "excludeNodes",
                                 doc="Allow the user to specify a set of nodes to be excluded from damage.")
    freezeDamage = PYB11property("bool", "freezeDamage", "freezeDamage",
                                 doc="Optionally freeze damage magnitude on each point")

#-------------------------------------------------------------------------------
# Add the restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, DamageModel)
