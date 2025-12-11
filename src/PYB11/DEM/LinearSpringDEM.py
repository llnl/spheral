#-------------------------------------------------------------------------------
# SPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from DEMBase import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralDEM")
class LinearSpringDEM(DEMBase):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using TimeStepType = typename DEMBase<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
"""
    
    def pyinit(self,
               dataBase = "const DataBase<%(Dimension)s>&",
               normalSpringConstant = "const Scalar",
               normalRestitutionCoefficient = "const Scalar",
               tangentialSpringConstant = "const Scalar",
               tangentialRestitutionCoefficient = "const Scalar",
               dynamicFrictionCoefficient = "const Scalar",
               staticFrictionCoefficient = "const Scalar",
               rollingFrictionCoefficient = "const Scalar",
               torsionalFrictionCoefficient = "const Scalar",
               cohesiveTensileStrength = "const Scalar",
               shapeFactor = "const Scalar",
               stepsPerCollision = "const Scalar",
               enableFastTimeStepping = "const bool",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "DEMBase constructor"

    @PYB11virtual
    @PYB11const
    def dt(self,
           dataBase = "const DataBase<%(Dimension)s>&", 
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

    @PYB11const
    def variableTimeStep(self,
                         dataBase = "const DataBase<%(Dimension)s>&", 
                         state = "const State<%(Dimension)s>&",
                         derivs = "const StateDerivatives<%(Dimension)s>&",
                         currentTime = "const Scalar"):
        "dt approach with fast time stepping enabled."
        return "TimeStepType"

    @PYB11const
    def fixedTimeStep(self):
        "simple dt based on contact time for when fast time stepping is off."
        return "TimeStepType"

    def recomputeContactDuration(self):
        "recalculates the contact time for use in the time step."
        return "void"

    @PYB11virtual
    def initializeProblemStartup(self,
                                 dataBase = "DataBase<%(Dimension)s>&"):
        "Tasks we do once on problem startup."
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
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "calculate the derivatives for Linear Spring DEM."
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

    @PYB11const
    def momentOfInertia(self,
                        massi = "const Scalar",
                        partialRadiusi = "const Scalar"):
        return "Scalar"


    def setMomentOfInertia(self):
        return "void"

    enableFastTimeStepping = PYB11property("bool", "enableFastTimeStepping", "enableFastTimeStepping", doc="activate fast time stepping")  
    normalSpringConstant = PYB11property("Scalar", "normalSpringConstant", "normalSpringConstant", doc="normal spring constant")
    normalRestitutionCoefficient = PYB11property("Scalar", "normalRestitutionCoefficient", "normalRestitutionCoefficient", doc="normal restitution coefficient")
    tangentialSpringConstant = PYB11property("Scalar", "tangentialSpringConstant", "tangentialSpringConstant", doc="tangential spring constant")
    tangentialRestitutionCoefficient = PYB11property("Scalar", "tangentialRestitutionCoefficient", "tangentialRestitutionCoefficient", doc="tangential restitution coefficient")
    cohesiveTensileStrength = PYB11property("Scalar", "cohesiveTensileStrength", "cohesiveTensileStrength", doc="constant for normal cohesion constant")
    
    dynamicFrictionCoefficient = PYB11property("Scalar", "dynamicFrictionCoefficient", "dynamicFrictionCoefficient", doc="sliding friction coefficient - dynamic")
    staticFrictionCoefficient = PYB11property("Scalar", "staticFrictionCoefficient", "staticFrictionCoefficient", doc="sliding friction coefficient - static")
    rollingFrictionCoefficient = PYB11property("Scalar", "rollingFrictionCoefficient", "rollingFrictionCoefficient", doc="rolling friction coefficient")
    torsionalFrictionCoefficient = PYB11property("Scalar", "torsionalFrictionCoefficient", "torsionalFrictionCoefficient", doc="torsional friction coefficient")
    
    shapeFactor = PYB11property("Scalar", "shapeFactor", "shapeFactor", doc="shape factor - simple approach to non-spherical particles")
    normalBeta = PYB11property("Scalar", "normalBeta", "normalBeta", doc="a damping parameter")
    tangentialBeta = PYB11property("Scalar", "tangentialBeta", "tangentialBeta", doc="a damping parameter")
    collisionDuration = PYB11property("Scalar", "collisionDuration", "collisionDuration", doc="duration of a contact")

    momentOfInertia = PYB11property("const FieldList<%(Dimension)s, Scalar>&","momentOfInertia", returnpolicy="reference_internal")
    maximumOverlap = PYB11property("const FieldList<%(Dimension)s, Scalar>&","maximumOverlap", returnpolicy="reference_internal")
    newMaximumOverlap = PYB11property("const FieldList<%(Dimension)s, Scalar>&","newMaximumOverlap", returnpolicy="reference_internal")
