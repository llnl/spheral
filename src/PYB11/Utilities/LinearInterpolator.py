#-------------------------------------------------------------------------------
# LinearInterpolator
#-------------------------------------------------------------------------------
from PYB11Generator import *

class LinearInterpolator:

    """Encapsulates the algorithm and data for linear regression in 1D
Assumes the results is interpolated as y_interp = a*x + b"""

    def pyinit(self):
        "Default constuctor -- returns a non-functional interpolator until initialized"
        return

    def pyinit_func(self,
                    xmin = "double",
                    xmax = "double",
                    n = "size_t",
                    F = "const PythonBoundFunctors::SpheralFunctor<double, double>&"):
        "Constructs an interpolator based on the given function sampled in x in [xmin, xmax]"
        return

    def pyinit_vals(self,
                    xvals = "const std::vector<double>&",
                    yvals = "const std::vector<double>&"):
        "Constructs an interpolator for yvals sampled in x in [xmin, xmax]"
        return

    def initialize(self,
                   xmin = "double",
                   xmax = "double",
                   n = "size_t",
                   F = "const PythonBoundFunctors::SpheralFunctor<double, double>&"):
        "Initializes the interpolator based on the given function sampled in x in [xmin, xmax]"
        return "void"

    @PYB11pycppname("initialize")
    def initialize_vals(self,
                        xvals = "const std::vector<double>&",
                        yvals = "const std::vector<double>&"):
        "Initializes the interpolator for xvals, yvals"
        return "void"

    @PYB11const
    def __call__(self,
                 x = "const double"):
        "Returns the interpolated value <y>(x)"
        return "double"

    # Attributes
    slope = PYB11property(doc="Fitted slope (a) for y = a*x + b")
    yintercept = PYB11property(doc="Fitted y-intercept (b) for y = a*x + b")
