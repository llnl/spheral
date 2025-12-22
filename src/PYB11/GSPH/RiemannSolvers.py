from PYB11Generator import *
#from RiemannSolverBaseAbstractMethods import *
#-------------------------------------------------------------------------------
# Base class for riemann solver wave speeds
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralGSPH")
class RiemannSolverBase:

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    """

    def pyinit(slopeLimiter = "LimiterBase<%(Dimension)s>&",
               waveSpeed = "WaveSpeedBase<%(Dimension)s>&",
               linearReconstruction = "const bool"):
        "slope limiter constructor"


    linearReconstruction = PYB11property("bool", "linearReconstruction", "linearReconstruction", doc="linear reconstruction of L and R state for riemann problem") 

    waveSpeed = PYB11property("WaveSpeedBase<%(Dimension)s>&", "waveSpeed",returnpolicy="reference_internal", doc="wave speed object")
    limiter = PYB11property("LimiterBase<%(Dimension)s>&", "limiter",returnpolicy="reference_internal", doc="slope limiter object")

#-------------------------------------------------------------------------------
# HLLC Approximate Riemann Solver
#-------------------------------------------------------------------------------
class HLLC(RiemannSolverBase):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    """

    def pyinit(slopeLimiter = "LimiterBase<%(Dimension)s>&",
               waveSpeed = "WaveSpeedBase<%(Dimension)s>&",
               linearReconstruction = "const bool"):
        "slope limiter constructor"


#-------------------------------------------------------------------------------
# HLLC Approximate Riemann Solver
#-------------------------------------------------------------------------------
class SecondOrderArtificialViscosity(RiemannSolverBase):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    """

    def pyinit(Cl = "const Scalar",
               Cq = "const Scalar",
               slopeLimiter = "LimiterBase<%(Dimension)s>&",
               waveSpeed = "WaveSpeedBase<%(Dimension)s>&",
               linearReconstruction = "const bool"):
        "slope limiter constructor"

    Cl = PYB11property("Scalar", "Cl", "Cl", doc="linear artificial viscosity coefficient") 
    Cq = PYB11property("Scalar", "Cq", "Cq", doc="quadratic artificial viscosity coefficient") 
