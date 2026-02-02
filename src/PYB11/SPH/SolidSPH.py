#-------------------------------------------------------------------------------
# SolidSPH
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SPHBase import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralSPH")
@PYB11dynamic_attr
class SolidSPH(SPHBase):
    "SolidSPH -- The SPH/ASPH solid material hydrodynamic package for Spheral++."

    PYB11typedefs = """
  using Scalar = typename %(Dimension)s::Scalar;
  using Vector = typename %(Dimension)s::Vector;
  using Tensor = typename %(Dimension)s::Tensor;
  using SymTensor = typename %(Dimension)s::SymTensor;
  using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
  using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
"""
    
    def pyinit(dataBase = "DataBase<%(Dimension)s>&",
               Q = "ArtificialViscosityHandle<%(Dimension)s>&",
               W = "const TableKernel<%(Dimension)s>&",
               WPi = "const TableKernel<%(Dimension)s>&",
               WGrad = "const TableKernel<%(Dimension)s>&",
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
               damageRelieveRubble = "const bool",
               strengthInDamage = "const bool",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "SolidSPH constructor"

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
    # Properties
    GradKernel = PYB11property("const TableKernel<%(Dimension)s>&", "GradKernel",
                               doc="Kernel for estimating velocity gradient")
    damageRelieveRubble = PYB11property("bool", "damageRelieveRubble", "damageRelieveRubble", 
                                        doc="Control whether allow damaged material to have stress relieved.")
    strengthInDamage = PYB11property("bool", "strengthInDamage", "strengthInDamage", 
                                     doc="Should we allow damaged material to support strength?")

    DdeviatoricStressDt =  PYB11property("const FieldList<%(Dimension)s, SymTensor>&", "DdeviatoricStressDt", returnpolicy="reference_internal")
    bulkModulus =          PYB11property("const FieldList<%(Dimension)s, Scalar>&",    "bulkModulus",         returnpolicy="reference_internal")
    shearModulus =         PYB11property("const FieldList<%(Dimension)s, Scalar>&",    "shearModulus",        returnpolicy="reference_internal")
    yieldStrength =        PYB11property("const FieldList<%(Dimension)s, Scalar>&",    "yieldStrength",       returnpolicy="reference_internal")
    plasticStrain0 =       PYB11property("const FieldList<%(Dimension)s, Scalar>&",    "plasticStrain0",      returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, SolidSPH)
