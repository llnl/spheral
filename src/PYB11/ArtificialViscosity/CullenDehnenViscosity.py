#-------------------------------------------------------------------------------
# CullenDehnenViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class CullenDehnenViscosity(Physics):
    """An implementation of the Cullen Dehnen Viscosity, with Hopkins alterations
Computes a correction for the scaling.
References:
Cullen, L., & Dehnen, W. 2010, MNRAS, 408, 669
Hopkins arXiv:1409.7395
"""

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ThirdRankTensor = typename %(Dimension)s::ThirdRankTensor;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               W = "const TableKernel<%(Dimension)s>&",
               alphMax = ("const Scalar", "2.0"),
               alphMin = ("const Scalar", "0.02"),
               betaC = ("const Scalar", "0.7"),
               betaD = ("const Scalar", "0.05"),
               betaE = ("const Scalar", "1.0"),
               fKern = ("const Scalar", "0.3333333333"),
               boolHopkins = ("const bool", "true")):
        "Cullen-Dehnen viscosity evolution constructor"

    #...........................................................................
    # Virtual methods
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
    def initializeProblemStartup(self,
                                 dataBase = "DataBase<%(Dimension)s>&"):
        "An optional hook to initialize once when the problem is starting up."
        return "void"

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
    def finalizeDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Provide a hook to be called after all physics packages have had their evaluateDerivatives method called, but before anyone does anything with those derivatives."
        return "void"

    #...........................................................................
    # Properties
    alphMax = PYB11property("Scalar", "alphMax", "alphMax")
    alphMin = PYB11property("Scalar", "alphMin", "alphMin")
    betaE = PYB11property("Scalar", "betaE", "betaE")
    betaD = PYB11property("Scalar", "betaD", "betaD")
    betaC = PYB11property("Scalar", "betaC", "betaC")
    fKern = PYB11property("Scalar", "fKern", "fKern")
    boolHopkins = PYB11property("bool", "boolHopkins", "boolHopkins")
    kernel = PYB11property("const TableKernel<%(Dimension)s>&", "kernel", returnpolicy="reference_internal")
    ClMultiplier = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "ClMultiplier",
                                 doc="Correction multiplier for the linear term")
    CqMultiplier = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "CqMultiplier",
                                 doc="Correction multiplier for the quadratic term")
    PrevDvDt = PYB11property("const FieldList<%(Dimension)s, Vector>&", "PrevDvDt", returnpolicy="reference_internal")
    PrevDivV = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "PrevDivV", returnpolicy="reference_internal")
    PrevDivV2 = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "PrevDivV2", returnpolicy="reference_internal")
    CullAlpha = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "CullAlpha", returnpolicy="reference_internal")
    CullAlpha2  = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "CullAlpha2", returnpolicy="reference_internal")
    DalphaDt = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "DalphaDt", returnpolicy="reference_internal")
    alphaLocal = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "alphaLocal", returnpolicy="reference_internal")
    R = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "R", returnpolicy="reference_internal")
    vsig = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "vsig", returnpolicy="reference_internal")
    
#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, CullenDehnenViscosity, pure_virtual=False, virtual=True)
PYB11inject(RestartMethods, CullenDehnenViscosity)
