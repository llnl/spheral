from PYB11Generator import *
from WaveSpeedBaseAbstractMethods import *

#-------------------------------------------------------------------------------
# Base class for riemann solver wave speeds
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralGSPH")
class WaveSpeedBase:

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "wave speed constructor"

PYB11inject(WaveSpeedBaseAbstractMethods, WaveSpeedBase, pure_virtual=True)


#-------------------------------------------------------------------------------
# acoustic wave speed
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralGSPH")
class AcousticWaveSpeed(WaveSpeedBase):

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "acoustic wave speed constructor"

    @PYB11virtual
    @PYB11const
    def waveSpeed(self,
                  rhoi = "const Scalar",
                  rhoj = "const Scalar",
                  ci   = "const Scalar",
                  cj   = "const Scalar",
                  ui   = "const Scalar",
                  uj   = "const Scalar",
                  Si   = "Scalar&",
                  Sj   = "Scalar&"):
        "calculate wave speed."
        return "void"

#-------------------------------------------------------------------------------
# Davis wave speed
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralGSPH")
class DavisWaveSpeed(WaveSpeedBase):

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "davis wave speed constructor"

    @PYB11virtual
    @PYB11const
    def waveSpeed(self,
                  rhoi = "const Scalar",
                  rhoj = "const Scalar",
                  ci   = "const Scalar",
                  cj   = "const Scalar",
                  ui   = "const Scalar",
                  uj   = "const Scalar",
                  Si   = "Scalar&",
                  Sj   = "Scalar&"):
        "calculate wave speed."
        return "void"


#-------------------------------------------------------------------------------
# einfeldt wave speed
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralGSPH")
class EinfeldtWaveSpeed(WaveSpeedBase):

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  """

    def pyinit():
        "einfeldt wave speed constructor"

    @PYB11virtual
    @PYB11const
    def waveSpeed(self,
                  rhoi = "const Scalar",
                  rhoj = "const Scalar",
                  ci   = "const Scalar",
                  cj   = "const Scalar",
                  ui   = "const Scalar",
                  uj   = "const Scalar",
                  Si   = "Scalar&",
                  Sj   = "Scalar&"):
        "calculate wave speed."
        return "void"

