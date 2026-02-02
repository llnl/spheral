#-------------------------------------------------------------------------------
# BiCubicInterpolator
#-------------------------------------------------------------------------------
from PYB11Generator import *
from XYInterpolator import *

class BiCubicInterpolator(XYInterpolator):

    """BiCubicInterpolator

Encapsulates the algorithm and data for cubic interpolation in 2D
Assumes the results is interpolated as
  <F(x,y)> = [1 x x^2 x^3][c00 c01 c02 c03][1]
                          [c10 c11 c12 c13][y]
                          [c20 c21 c22 c23][y^2]
                          [c30 c31 c32 c33][y^3]
Assumes we provide functors for F and gradF for fitting, where
     F(x,y) -> double
 gradF(x,y) -> Dim<2>::SymTensor

See https://en.wikipedia.org/wiki/Bicubic_interpolation"""

    def pyinit(self):
        "Default constuctor -- returns a non-functional interpolator until initialized"
        return

    def pyinit1(self,
                xmin = "const double",
                xmax = "const double",
                ymin = "const double",
                ymax = "const double",
                nx = "const size_t",
                ny = "const size_t",
                F = "const Spheral::PythonBoundFunctors::Spheral2ArgFunctor<double, double, double>&",
                xlog = ("const bool", "false"),
                ylog = ("const bool", "false")):
        """Returns an interpolator for z-vals sampled in (x,y) in ([xmin, xmax], [ymin, ymax]).
In this case without a gradient function, we sample a stencil of 16 values
in each cell for fitting."""
        return

    def pyinit2(self,
                xmin = "const double",
                xmax = "const double",
                ymin = "const double",
                ymax = "const double",
                nx = "const size_t",
                ny = "const size_t",
                F = "const Spheral::PythonBoundFunctors::Spheral2ArgFunctor<double, double, double>&",
                Fgrad = "const Spheral::PythonBoundFunctors::Spheral2ArgFunctor<double, double, Dim<2>::SymTensor>&",
                xlog = ("const bool", "false"),
                ylog = ("const bool", "false")):
        """Returns an interpolator for z-vals sampled in (x,y) in ([xmin, xmax], [ymin, ymax]).
Since this form has the gradient, we only need to sample the function at the
four corners of each interpolation cell."""
        return

    @PYB11const
    def __call__(self,
                 x = "const double",
                 y = "const double"):
        "Returns the interpolated value <F>(x,y)"
        return "double"

    @PYB11const
    def prime_x(self,
                x = "const double",
                y = "const double"):
        """Returns the interpolated value partial_x <F>(x,y)"""
        return "double"

    @PYB11const
    def prime_y(self,
                x = "const double",
                y = "const double"):
        """Returns the interpolated value partial_y <F>(x,y)"""
        return "double"

    @PYB11const
    def prime2_xx(self,
                  x = "const double",
                  y = "const double"):
        """Returns the interpolated value partial_xx <F>(x,y)"""
        return "double"

    @PYB11const
    def prime2_xy(self,
                  x = "const double",
                  y = "const double"):
        """Returns the interpolated value partial_xy <F>(x,y)"""
        return "double"

    @PYB11const
    def prime2_yy(self,
                  x = "const double",
                  y = "const double"):
        """Returns the interpolated value partial_yy <F>(x,y)"""
        return "double"
