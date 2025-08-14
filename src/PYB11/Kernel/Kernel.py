#-------------------------------------------------------------------------------
# Generic Kernel bindings.
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension", "Descendant")
class Kernel:

    PYB11typedefs = """
    using Vector = typename %(Dimension)s::Vector;
    using SymTensor = typename %(Dimension)s::SymTensor;
"""

    @PYB11pyname("__call__")
    @PYB11cppname("operator()")
    @PYB11const
    def __call__etaij(self,
                      etaij = "const double&",
                      Hdet = "const double&"):
        "Return the kernel value at the given normalized radius, where etaij = magnitude(H*(posi - posj)), Hdet = ||H||"
        return "double"

    @PYB11pyname("__call__")
    @PYB11cppname("operator()")
    @PYB11const
    def __call__etaj_etai(self,
                          etaj = "const Vector&",
                          etai = "const Vector&",
                          Hdet = "const double&"):
        "Return the kernel value for points at etaj and etai, where etaj = H*posj, etai = H*posi, Hdet = ||H||"
        return "double"

    #...........................................................................
    @PYB11const
    def grad(self,
             etaij = "const double&",
             Hdet = "const double&"):
        """Return the gradient value (scalar) at the given normalized radius, where etaij = magnitude(H*(posi - posj)), Hdet = ||H||.
The full gradient can be constructed from the result as \vec{grad} = H * \vec{etaij}/|\vec{etaij}| * result"""
        return "double"

    @PYB11pycppname("grad")
    @PYB11const
    def another_grad(self,
                     etaj = "const Vector&",
                     etai = "const Vector&",
                     Hdet = "const double&"):
        """Return the gradient (scalar) value for points at etaj and etai, where etaj = H*posj, etai = H*posi, Hdet = ||H||.
The full gradient can be constructed from the result as \vec{grad} = H * \vec{etaij}/|\vec{etaij}| * result"""
        return "double"

    #...........................................................................
    @PYB11const
    def grad2(self,
              etaij = "const double&",
              Hdet = "const double&"):
        "Return the second derivative value (scalar) at the given normalized radius, where etaij = magnitude(H*(posi - posj)), Hdet = ||H||"
        return "double"

    @PYB11pycppname("grad")
    @PYB11const
    def another_grad2(self,
                      etaj = "const Vector&",
                      etai = "const Vector&",
                      Hdet = "const double&"):
        """Return the gradient (scalar) value for points at etaj and etai, where etaj = H*posj, etai = H*posi, Hdet = ||H||.
The full gradient can be constructed from the result as \vec{grad} = H * \vec{etaij}/|\vec{etaij}| * result"""
        return "double"

    #...........................................................................
    @PYB11const
    def kernelValue(self,
                    etaij = "double",
                    Hdet = "const double"):
        "Return the value of the kernel"
        return "double"

    @PYB11const
    def gradValue(self,
                  etaij = "double",
                  Hdet = "const double"):
        "Return the value of the gradient of the kernel"
        return "double"

    @PYB11const
    def grad2Value(self,
                   etaij = "double",
                   Hdet = "const double"):
        "Return the value of the second derivative of the kernel"
        return "double"

    #...........................................................................
    # Protected methods
    @PYB11protected
    @PYB11const
    def setVolumeNormalization(self,
                               val = "double"):
        "Set the volume normalization constant."
        return "void"

    @PYB11protected
    @PYB11const
    def setKernelExtent(self,
                        val = "double"):
        "Set the volume normalization constant."
        return "void"

    @PYB11protected
    @PYB11const
    def setInflectionPoint(self,
                           val = "double"):
        "Set the volume normalization constant."
        return "void"

    #...........................................................................
    # Properties
    volumeNormalization = PYB11property("double", "volumeNormalization", doc="The volume normalization to ensure the integral of the kennel is unity.")
    kernelExtent = PYB11property("double", "kernelExtent", doc="The spatial radius of the kernel in dimensionless (x/h) units.")
    inflectionPoint = PYB11property("double", "inflectionPoint", doc="The inflection point radius in dimensionless (x/h) units.")

#-------------------------------------------------------------------------------
# BSpline
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class BSplineKernel(Kernel):
    PYB11typedefs=""

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# W4Spline
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class W4SplineKernel(Kernel):
    PYB11typedefs=""

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# Gaussian
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class GaussianKernel(Kernel):
    PYB11typedefs=""

    def pyinit(self,
               extent = "double"):
        "Constructor with the given eta cutoff"

#-------------------------------------------------------------------------------
# SuperGaussian
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class SuperGaussianKernel(Kernel):
    PYB11typedefs=""

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# PiGaussianSpline
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class PiGaussianKernel(Kernel):
    PYB11typedefs=""

    def pyinit(self):
        "Default constructor"

    def pyinit(self,
               K = "double"):
        "Construct with K"

    K = PYB11property("double", "getK", "setK", "Set the K constant")

#-------------------------------------------------------------------------------
# Hat
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class HatKernel(Kernel):
    PYB11typedefs=""

    def pyinit(self,
               eta0 = "double",
               W0 = "double"):
        "Construct with (eta0, K)"
    
    eta0 = PYB11property("double", "eta0", doc="eta0")
    W0 = PYB11property("double", "W0", doc="W0")

#-------------------------------------------------------------------------------
# Sinc
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class SincKernel(Kernel):
    PYB11typedefs=""

    def pyinit(self, extent="double"):
        "Construct with the given extent in eta"

#-------------------------------------------------------------------------------
# NSincPolynomial
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class NSincPolynomialKernel(Kernel):
    PYB11typedefs=""

    def pyinit(self, order="int"):
        "Construct with the specified sinc order"

#-------------------------------------------------------------------------------
# NBSpline
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class NBSplineKernel(Kernel):
    PYB11typedefs=""

    def pyinit(self, order="int"):
        "Construct using the given order of b-spline"

    @PYB11const
    def factorial(self, n="int"):
        "Factorial n!"
        return "int"

    @PYB11const
    def binomialCoefficient(self, n="int", m="int"):
        "Compute the binomial coefficient of (n,m)"
        return "int"

    @PYB11const
    def oneSidedPowerFunction(self, s="double", exponent="int"):
        return "double"

    order = PYB11property("int", "order", "setOrder", "The order of the kernel")

#-------------------------------------------------------------------------------
# Table
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11template_dict({"Descendant": "TableKernelView<%(Dimension)s>"})
class TableKernelView(Kernel):
    PYB11typedefs=""

    def pyinit(self):
        "Default constructor"

@PYB11template("Dimension")
class TableKernel(TableKernelView):
    PYB11typedefs="""
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using SymTensor = typename %(Dimension)s::SymTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               kernel = "const BSplineKernel<%(Dimension)s>&",
               numPoints = ("const unsigned", "100"),
               minNperh = ("const double", "0.25"),
               maxNperh = ("const double", "64.0")):
        "Construct with BSpline kernel"

    def pyinita(self,
                kernel = "const W4SplineKernel<%(Dimension)s>&",
                numPoints = ("const unsigned", "100"),
                minNperh = ("const double", "0.25"),
                maxNperh = ("const double", "64.0")):
        "Construct with W4Spline kernel"

    def pyinitb(self,
                kernel = "const GaussianKernel<%(Dimension)s>&",
                numPoints = ("const unsigned", "100"),
                minNperh = ("const double", "0.25"),
                maxNperh = ("const double", "64.0")):
        "Construct with Gaussian kernel"

    def pyinitc(self,
                kernel = "const SuperGaussianKernel<%(Dimension)s>&",
                numPoints = ("const unsigned", "100"),
                minNperh = ("const double", "0.25"),
                maxNperh = ("const double", "64.0")):
        "Construct with SuperGaussian kernel"

    def pyinitd(self,
                kernel = "const PiGaussianKernel<%(Dimension)s>&",
                numPoints = ("const unsigned", "100"),
                minNperh = ("const double", "0.25"),
                maxNperh = ("const double", "64.0")):
        "Construct with PiGaussian kernel"

    def pyinite(self,
                kernel = "const HatKernel<%(Dimension)s>&",
                numPoints = ("const unsigned", "100"),
                minNperh = ("const double", "0.25"),
                maxNperh = ("const double", "64.0")):
        "Construct with Hat kernel"

    def pyinitf(self,
                kernel = "const SincKernel<%(Dimension)s>&",
                numPoints = ("const unsigned", "100"),
                minNperh = ("const double", "0.25"),
                maxNperh = ("const double", "64.0")):
        "Construct with Sinc kernel"

    def pyinitg(self,
                kernel = "const NSincPolynomialKernel<%(Dimension)s>&",
                numPoints = ("const unsigned", "100"),
                minNperh = ("const double", "0.25"),
                maxNperh = ("const double", "64.0")):
        "Construct with NSincPolynomial kernel"

    def pyinith(self,
                kernel = "const QuarticSplineKernel<%(Dimension)s>&",
                numPoints = ("const unsigned", "100"),
                minNperh = ("const double", "0.25"),
                maxNperh = ("const double", "64.0")):
        "Construct with Quartic spline kernel"

    def pyiniti(self,
                kernel = "const QuinticSplineKernel<%(Dimension)s>&",
                numPoints = ("const unsigned", "100"),
                minNperh = ("const double", "0.25"),
                maxNperh = ("const double", "64.0")):
        "Construct with Quintic spline kernel"

    def pyinitj(self,
                kernel = "const NBSplineKernel<%(Dimension)s>&",
                numPoints = ("const unsigned", "100"),
                minNperh = ("const double", "0.25"),
                maxNperh = ("const double", "64.0")):
        "Construct with NBSpline kernel"

    def pyinitk(self,
                kernel = "const WendlandC2Kernel<%(Dimension)s>&",
                numPoints = ("const unsigned", "100"),
                minNperh = ("const double", "0.25"),
                maxNperh = ("const double", "64.0")):
        "Construct with WendlandC2 kernel"

    def pyinitl(self,
                kernel = "const WendlandC4Kernel<%(Dimension)s>&",
                numPoints = ("const unsigned", "100"),
                minNperh = ("const double", "0.25"),
                maxNperh = ("const double", "64.0")):
        "Construct with WendlandC4 kernel"

    def pyinitm(self,
                kernel = "const WendlandC6Kernel<%(Dimension)s>&",
                numPoints = ("const unsigned", "100"),
                minNperh = ("const double", "0.25"),
                maxNperh = ("const double", "64.0")):
        "Construct with WendlandC6 kernel"

    #...........................................................................
    # Methods
    @PYB11const
    @PYB11implementation("""[](const TableKernel<%(Dimension)s>& self, const Vector& etaj, const Vector& etai, const SymTensor& H) -> py::tuple {
        Scalar W, deltaWsum;
        Vector gradW;
        self.kernelAndGrad(etaj, etai, H, W, gradW, deltaWsum);
        return py::make_tuple(W, gradW, deltaWsum);
      }""")
    def kernelAndGrad(self,
                      etaj = "const Vector&",
                      etai = "const Vector&",
                      H = "const SymTensor&"):
        "Simultaneously compute the W(etaj, etai), gradW(etaj, etai), deltaWsum(etaj, etai) -- returns a tuple of those results."
        return "py::tuple"

    @PYB11const
    @PYB11implementation("""[](const TableKernel<%(Dimension)s>& self, const Scalar& etaij, const Scalar& Hdet) -> py::tuple {
        Scalar W, gW;
        self.kernelAndGradValue(etaij, Hdet, W, gW);
        return py::make_tuple(W, gW);
      }""")
    def kernelAndGradValue(self,
                           etaij = "const Scalar",
                           Hdet = "const Scalar"):
        return "py::tuple"

    @PYB11const
    def kernelAndGradValues(self,
                            etaij = "const std::vector<Scalar>&",
                            Hdets = "const std::vector<Scalar>&",
                            kernelValues = "std::vector<Scalar>&",
                            gradValues = "std::vector<Scalar>&"):
        return "void"

    @PYB11const
    def kernelValueSPH(self,
                       etaij = "const Scalar"):
        "Compute the kernel value appropriate for use in the SPH variable h 'ideal h' calculation"
        return "Scalar"

    @PYB11const
    def kernelValueASPH(self,
                        etaij = "const Scalar",
                        nPerh = "const Scalar"):
        "Compute the kernel value appropriate for use in the ASPH variable H 'ideal H' calculation"
        return "Scalar"

    @PYB11const
    def equivalentNodesPerSmoothingScale(self,
                                         Wsum = "Scalar"):
        "Compute the nPerh that corresponds to the Wsum value"
        return "Scalar"

    @PYB11const
    def equivalentWsum(self,
                       nPerh = "Scalar"):
        "Compute the Wsum that corresponds to the  nPerh value"
        return "Scalar"

    #...........................................................................
    # Properties
    numPoints = PYB11property("size_t", doc="The number of points in the table")
    minNperhLookup = PYB11property("double", doc="The lower limit for looking up the effective nPerh")
    maxNperhLookup = PYB11property("double", doc="The upper limit for looking up the effective nPerh")
    Winterpolator = PYB11property(doc = "W(x) interpolator")
    gradWinterpolator = PYB11property(doc = "grad W(x) interpolator")
    grad2Winterpolator = PYB11property(doc = "grad^2 W(x) interpolator")
    nPerhInterpolator = PYB11property(doc = "nperh(x) interpolator (SPH)")
    WsumInterpolator = PYB11property(doc = "Wsum(x) interpolator (SPH)")

#-------------------------------------------------------------------------------
# WendlandC2
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class WendlandC2Kernel(Kernel):
    PYB11typedefs=""

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# WendlandC4
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class WendlandC4Kernel(Kernel):
    PYB11typedefs=""

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# WendlandC6
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class WendlandC6Kernel(Kernel):
    PYB11typedefs=""

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# QuarticSpline
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class QuarticSplineKernel(Kernel):
    PYB11typedefs=""

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# QuinticSpline
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class QuinticSplineKernel(Kernel):
    PYB11typedefs=""

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# ExpInv
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class ExpInvKernel(Kernel):
    PYB11typedefs=""

    def pyinit(self):
        "Default constructor"

