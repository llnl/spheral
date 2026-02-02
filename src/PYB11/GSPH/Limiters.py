from PYB11Generator import *
from LimiterBaseAbstractMethods import *

#-------------------------------------------------------------------------------
# Base class for riemann solver wave speeds
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralGSPH")
class LimiterBase:

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit(TVD = "bool",
               symmetric = "bool"):
        "slope limiter constructor"

    @PYB11virtual
    @PYB11const
    def slopeLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def isTVD(self):
        "is this a TVD limiter?."
        return "bool"

    @PYB11virtual
    @PYB11const
    def isSymmetric(self):
        "is this a Symmetric limiter?."
        return "bool"

PYB11inject(LimiterBaseAbstractMethods, LimiterBase, pure_virtual=True)

#-------------------------------------------------------------------------------
# minmod limiter
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralGSPH")
class MinModLimiter(LimiterBase):

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "minmod slope limiter constructor"

    @PYB11virtual
    @PYB11const
    def slopeLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def fluxLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

#-------------------------------------------------------------------------------
# VanLeer limiter
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralGSPH")
class VanLeerLimiter(LimiterBase):

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "VanLeer slope limiter constructor"

    @PYB11virtual
    @PYB11const
    def slopeLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def fluxLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

#-------------------------------------------------------------------------------
# VanAlba limiter
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralGSPH")
class VanAlbaLimiter(LimiterBase):

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "VanAlba slope limiter constructor"

    @PYB11virtual
    @PYB11const
    def slopeLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def fluxLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

#-------------------------------------------------------------------------------
# Superbee limiter
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralGSPH")
class SuperbeeLimiter(LimiterBase):

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "Superbee slope limiter constructor"

    @PYB11virtual
    @PYB11const
    def slopeLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def fluxLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

#-------------------------------------------------------------------------------
# Ospre limiter
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralGSPH")
class OspreLimiter(LimiterBase):

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "Ospre slope limiter constructor"

    @PYB11virtual
    @PYB11const
    def slopeLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def fluxLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

#-------------------------------------------------------------------------------
# Barth Jespersen limiter
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralGSPH")
class BarthJespersenLimiter(LimiterBase):

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "Barth Jespersen slope limiter constructor"

    @PYB11virtual
    @PYB11const
    def slopeLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"

    @PYB11virtual
    @PYB11const
    def fluxLimiter(self,
                     x = "const Scalar"):
        "slope limiter from flux limiter."
        return "Scalar"


