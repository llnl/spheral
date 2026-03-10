#-------------------------------------------------------------------------------
# CubicHermiteInterpolator
#-------------------------------------------------------------------------------
from PYB11Generator import *

class CubicHermiteInterpolator:

    "An (optionally monotonic) form of cubic Hermite interpolation."

    def pyinit(self):
        "Default constructor -- returns a non-functional interpolator until initialized"
        return

    def pyinit_func(self,
                    xmin = "const double",
                    xmax = "const double",
                    n = "const size_t",
                    F = "const PythonBoundFunctors::SpheralFunctor<double, double>&"):
        "Constructs an interpolator based on the given function"
        return

    def pyinit_gradfunc(self,
                        xmin = "const double",
                        xmax = "const double",
                        n = "const size_t",
                        F = "const PythonBoundFunctors::SpheralFunctor<double, double>&",
                        Fgrad = "const PythonBoundFunctors::SpheralFunctor<double, double>&"):
        "Constructs an interpolator based on the given function and its gradient"
        return

    @PYB11implementation("[](const double xmin, const double xmax, py::list& yvals) { std::vector<double> yvec; for (auto y: yvals) yvec.push_back(y.cast<double>()); return new CubicHermiteInterpolator(xmin, xmax, yvec); }")
    def pyinit_table(self,
                     xmin = "const double",
                     xmax = "const double",
                     yvals = "py::list"):
        "Initialize from tabulated values"
        return

    @PYB11pycppname("initialize")
    def initialize_func(self,
                        xmin = "const double",
                        xmax = "const double",
                        n = "const size_t",
                        F = "const PythonBoundFunctors::SpheralFunctor<double, double>&"):
        "Initialize based on the given function"
        return "void"

    @PYB11pycppname("initialize")
    def initialize_gradfunc(self,
                            xmin = "const double",
                            xmax = "const double",
                            n = "const size_t",
                            F = "const PythonBoundFunctors::SpheralFunctor<double, double>&",
                            Fgrad = "const PythonBoundFunctors::SpheralFunctor<double, double>&"):
        "Initialize based on the given function and its gradient"
        return "void"

    @PYB11pycppname("initialize")
    def initialize_table(self,
                         xmin = "const double",
                         xmax = "const double",
                         yvals = "const std::vector<double>&"):
        "Initializes the interpolator for yvals sampled in x in [xmin, xmax]"
        return "void"

    def makeMonotonic(self):
        """Force interpolation to be monotonic.  This generally degrades accuracy, and can introduce structure between
tabulated knots (although that structure should still be monotonic between the tabulated function values)."""
        return "void"

    @PYB11const
    def __call__(self,
                 x = "const double"):
        "Returns the interpolated value <y>(x)"
        return "double"

    @PYB11pyname("__call__")
    @PYB11const
    def __call__i0(self,
                   x = "const double",
                   i0 = "const size_t"):
        "Returns the interpolated value <y>(x)"
        return "double"

    @PYB11cppname("operator[]")
    def __getitem__(self, index="size_t"):
        return "double"

    @PYB11const
    def prime(self,
              x = "const double"):
        "Return the interpolated derivative <dy/dx>(x)"
        return "double"

    @PYB11pycppname("prime")
    @PYB11const
    def prime_i0(self,
                 x = "const double",
                 i0 = "const size_t"):
        "Return the interpolated derivative <dy/dx>(x)"
        return "double"

    @PYB11const
    def prime2(self,
               x = "const double"):
        "Return the interpolated second derivative <d^2y/dx^2>(x)"
        return "double"

    @PYB11pycppname("prime2")
    @PYB11const
    def prime2_i0(self,
                  x = "const double",
                  i0 = "const size_t"):
        "Return the interpolated second derivative <d^2y/dx^2>(x)"
        return "double"

    @PYB11const
    def lowerBound(self,
                   x = "const double"):
        "Return the lower bound index in the table for the given x coordinate"
        return "size_t"

    @PYB11const
    def h00(self,
            t = "const double"):
        "h00 Hermite basis function"
        return "double"

    @PYB11const
    def h10(self,
            t = "const double"):
        "h10 Hermite basis function"
        return "double"

    @PYB11const
    def h01(self,
            t = "const double"):
        "h01 Hermite basis function"
        return "double"

    @PYB11const
    def h11(self,
            t = "const double"):
        "h11 Hermite basis function"
        return "double"

    # Attributes
    size = PYB11property(doc="The number of the tabulated values used")
    xmin = PYB11property(doc="Minimum x coordinate for table")
    xmax = PYB11property(doc="Maximum x coordinate for table")
    xstep = PYB11property(doc="delta x between tabulated values")

