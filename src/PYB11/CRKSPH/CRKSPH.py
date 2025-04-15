#-------------------------------------------------------------------------------
# CRKSPHBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from CRKSPHBase import *

@PYB11template("Dimension")
@PYB11module("SpheralCRKSPH")
@PYB11dynamic_attr
class CRKSPH(CRKSPHBase):
    "CRKSPHBase -- The CRKSPH/ACRKSPH hydrodynamic package for Spheral++."

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ThirdRankTensor = typename %(Dimension)s::ThirdRankTensor;
    using FourthRankTensor = typename %(Dimension)s::FourthRankTensor;
    using FifthRankTensor = typename %(Dimension)s::FifthRankTensor;
    using FacetedVolume = typename %(Dimension)s::FacetedVolume;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
    using PairAccelerationsType = typename CRKSPH<%(Dimension)s>::PairAccelerationsType;
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
               nTensile = "const double"):
        "Constructor"

    #...........................................................................
    # Virtual methods
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

    #...........................................................................
    # Properties
    pairAccelerations = PYB11property("const PairAccelerationsType&", "pairAccelerations", returnpolicy="reference_internal")
