#-------------------------------------------------------------------------------
# wave Speed pure virtual interface
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11ignore
class GenericRiemannHydroAbstractMethods:

    @PYB11virtual
    @PYB11const
    def firstDerivativesLoop(time = "const Scalar",
                             dt = "const Scalar",
                             dataBase = "const DataBase<%(Dimension)s>&",
                             state = "const State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        """main loop, evaluates the derivatives for the principle hydro mass density, velocity, and specific thermal energy."""
        return "void"
    
    @PYB11virtual
    @PYB11const
    def secondDerivativesLoop(time = "const Scalar",
                              dt = "const Scalar",
                              dataBase = "const DataBase<%(Dimension)s>&",
                              state = "const State<%(Dimension)s>&",
                              derivs = "StateDerivatives<%(Dimension)s>&"):
        """used to preprocess gradients need for main loop."""
        return "void"
