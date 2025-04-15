#-------------------------------------------------------------------------------
# JohnsonCookDamage
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from RestartMethods import *

@PYB11template("Dimension")
class JohnsonCookDamage(Physics):
    "JohnsonCookDamage -- an implementation of a Johnson-Cook damage law."

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
               D1 = "const Field<%(Dimension)s, Scalar>&",
               D2 = "const Field<%(Dimension)s, Scalar>&",
               D3 = "const double",
               D4 = "const double",
               D5 = "const double",
               epsilondot0 = "const double",
               Tcrit = "const double",
               sigmamax = "const double",
               efailmin = "const double"):
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

    #...........................................................................
    # Properties
    nodeList = PYB11property("const SolidNodeList<%(Dimension)s>&", returnpolicy="reference_internal")
    failureStrain = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    meltSpecificEnergy = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    D1 = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    D2 = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    D3 = PYB11property("double")
    D4 = PYB11property("double")
    D5 = PYB11property("double")
    epsilondot0 = PYB11property("double")
    Tcrit = PYB11property("double")
    sigmamax = PYB11property("double")
    efailmin = PYB11property("double")

#-------------------------------------------------------------------------------
# Add the restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, JohnsonCookDamage)
