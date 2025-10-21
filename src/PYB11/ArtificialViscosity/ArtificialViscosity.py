#-------------------------------------------------------------------------------
# ArtificialViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralArtificialViscosity")
class ArtificialViscosity(Physics):

    PYB11typedefs = """
  using Scalar = typename %(Dimension)s::Scalar;
  using Vector = typename %(Dimension)s::Vector;
  using Tensor = typename %(Dimension)s::Tensor;
  using SymTensor = typename %(Dimension)s::SymTensor;
  using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
  using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               Clinear = "const Scalar",
               Cquadratic = "const Scalar",
               kernel = "const TableKernel<%(Dimension)s>&"):
        "ArtificialViscosity constructor"

    #...........................................................................
    # Methods
    @PYB11pure_virtual
    @PYB11const
    def QPiTypeIndex(self):
        "Require ArtificialViscosities to specify the type_index of the descendant QPiType"
        return "std::type_index"

    @PYB11virtual
    @PYB11const
    def requireVelocityGradient(self):
        "Some AVs need the velocity gradient computed, so they should override this to true"
        return "bool"

    @PYB11virtual
    def updateVelocityGradient(self,
                               dataBase = "const DataBase<%(Dimension)s>&",
                               state = "const State<%(Dimension)s>&",
                               derivs = "const StateDerivatives<%(Dimension)s>&"):
        "Update the locally stored velocity gradient"
        return "void"

    @PYB11virtual
    def applyGhostBoundaries(self,
                             state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields"
        return "void"

    @PYB11virtual
    def initializeProblemStartup(self,
                                 dataBase = "DataBase<%(Dimension)s>&"):
        "Initialize the artificial viscosity for all FluidNodeLists in the given DataBase"
        return "void"

    @PYB11virtual
    def initializeProblemStartupDependencies(self,
                                             dataBase = "DataBase<%(Dimension)s>&",
                                             state = "State<%(Dimension)s>&",
                                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Initialize the artificial viscosity for all FluidNodeLists in the given DataBase"
        return "void"

    @PYB11virtual
    def initialize(self,
                   time = "const Scalar",
                   dt = "const Scalar",
                   dataBase = "const DataBase<%(Dimension)s>&",
                   state = "State<%(Dimension)s>&",
                   derivs = "StateDerivatives<%(Dimension)s>&"):
        "initialize is our chance to update the velocity gradient if needed"
        return "bool"

    @PYB11virtual
    def postStateUpdate(self,
                        time = "const Scalar",
                        dt = "const Scalar",
                        dataBase = "const DataBase<%(Dimension)s>&",
                        state = "State<%(Dimension)s>&",
                        derivs = "StateDerivatives<%(Dimension)s>&"):
        "Post-state update is our chance to update the velocity gradient if needed"
        return "bool"

    @PYB11const
    def curlVelocityMagnitude(self, DvDx="const Tensor&"):
        "Calculate the curl of the velocity given the stress tensor."
        return "Scalar"

    @PYB11const
    def calcBalsaraShearCorrection(self,
                                   DvDx = "const Tensor&",
                                   H = "const SymTensor&",
                                   cs = "const Scalar&"):
        "Find the Balsara shear correction multiplier"
        return "Scalar"

    @PYB11pure_virtual
    def updateManagedPtr(self):
        "Update member data for managed pointer."
        return "void"

    #...........................................................................
    # Properties
    Cl = PYB11property("Scalar", "Cl", "Cl",
                       doc="The linear coefficient")
    Cq = PYB11property("Scalar", "Cq", "Cq",
                       doc="The quadratic coefficient")
    balsaraShearCorrection = PYB11property("bool", "balsaraShearCorrection", "balsaraShearCorrection",
                                           doc="Toggle whether to use the Balsara suppression for shear flows")
    limiter = PYB11property("bool", "limiter", "limiter", doc="Defunct option, to be removed")
    epsilon2 = PYB11property("Scalar", "epsilon2", "epsilon2",
                             doc="Safety factor in denominator for Q")
    negligibleSoundSpeed = PYB11property("Scalar", "negligibleSoundSpeed", "negligibleSoundSpeed",
                                         doc="The negligible sound speed parameter for use in the limiter")
    maxViscousPressure = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "maxViscousPressure",
                                       doc="Store the maximum viscous pressure (Q) on each point")
    effViscousPressure = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "effViscousPressure",
                                       doc="Store the effective (averaged) viscous pressure (Q) on each point")
    DvDx = PYB11property("const FieldList<%(Dimension)s, Tensor>&", "DvDx",
                         doc="The velocity gradient used for AV")
    rigorousVelocityGradient = PYB11property("bool", "rigorousVelocityGradient", "rigorousVelocityGradient",
                                 doc="The multiplier for sound speed in the limiter")
    kernel = PYB11property("const TableKernel<%(Dimension)s>&", "kernel",
                           doc="Interpolation kernel for estimating velocty gradient (if needed)")
    
#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, ArtificialViscosity, pure_virtual=False)
PYB11inject(RestartMethods, ArtificialViscosity)
