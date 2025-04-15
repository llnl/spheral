#-------------------------------------------------------------------------------
# IvanoviSALEDamageModel
#-------------------------------------------------------------------------------
from PYB11Generator import *
from DamageModel import *
from RestartMethods import *

@PYB11template("Dimension")
class IvanoviSALEDamageModel(DamageModel):
    """The Ivanov damage model, hopefully close to how it's implemented in iSALE.
This damage model is most appropriate for rocky materials.

Refs:

Collins, G. S., Melosh, H. J., & Ivanov, B. A. (2004). Modeling damage and deformation in impact simulations.
  Meteoritics & Planetary Science, 39(2), 217. http://doi.wiley.com/10.1111/j.1945-5100.2004.tb00337.x

Raducan, S. D., Davison, T. M., Luther, R., & Collins, G. S. (2019). The role of asteroid strength, porosity and
  internal friction in impact momentum transfer. Icarus. https://doi.org/10.1016/J.ICARUS.2019.03.040

Lundborg, N. (1967). The strength-size relation of granite. International Journal of Rock Mechanics and
  Mining Sciences & Geomechanics Abstracts, 4(3):269
"""

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
               minPlasticFailure = "const double",
               plasticFailurePressureSlope = "const double",
               plasticFailurePressureOffset = "const double",
               tensileFailureStress = "const double",
               crackGrowthMultiplier = "const double",
               damageCouplingAlgorithm  = "const DamageCouplingAlgorithm",
               criticalDamageThreshold = "const double",
               mask = "const Field<%(Dimension)s, int>&"):
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
    # Properties
    minPlasticFailure = PYB11property("double", "minPlasticFailure",
                                      doc="Minimum plastic strain for shear failure")
    plasticFailurePressureSlope = PYB11property("double", "plasticFailurePressureSlope",
                                                doc="Slope of the failure plastic strain linear law")
    plasticFailurePressureOffset = PYB11property("double", "plasticFailurePressureOffset",
                                                 doc="Offset of the failure plastic strain linear law")
    tensileFailureStress = PYB11property("double", "tensileFailureStress",
                                         doc="The stress required for tensile failure")
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
PYB11inject(RestartMethods, IvanoviSALEDamageModel)
