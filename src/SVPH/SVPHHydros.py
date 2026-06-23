from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic SVPHHydro pattern.
#-------------------------------------------------------------------------------
FactoryString = """
class %(classname)s%(dim)s(%(basename)sBase%(dim)s):

    def __init__(self,
                 W,
                 Q,
                 cfl = 0.5,
                 useVelocityMagnitudeForDt = False,
                 useNewAccelerationMagnitudeForDt = False,
                 compatibleEnergyEvolution = True,
                 XSVPH = True,
                 linearConsistent = False,
                 generateVoid = True,
                 densityUpdate = RigorousSumDensity,
                 HUpdate = IdealH,
                 fcentroidal = 0.0,
                 fcellPressure = 0.0,
                 xmin = Vector%(dim)s(-1e100, -1e100, -1e100),
                 xmax = Vector%(dim)s( 1e100,  1e100,  1e100)):
        self._smoothingScaleMethod = %(smoothingScaleMethod)s%(dim)s()
        if xmin is None:
            xmin = Vector%(dim)s.one
        %(basename)sBase%(dim)s.__init__(self,
                                      self._smoothingScaleMethod,
                                      W,
                                      Q,
                                      cfl,
                                      useVelocityMagnitudeForDt,
                                      useNewAccelerationMagnitudeForDt,
                                      compatibleEnergyEvolution,
                                      XSVPH,
                                      linearConsistent,
                                      generateVoid,
                                      densityUpdate,
                                      HUpdate,
                                      fcentroidal,
                                      fcellPressure,
                                      xmin,
                                      xmax)
        return
"""

#-------------------------------------------------------------------------------
# Make 'em.
#-------------------------------------------------------------------------------
for ndim in dims:
    dim = "%id" % ndim
    for basename in ("SVPHFacetedHydro",):
        for prefix, smoothingScaleMethod in (("", "SPHSmoothingScale"),
                                             ("A", "ASPHSmoothingScale")):
            exec(FactoryString % {"dim"                  : dim,
                                  "basename"             : basename,
                                  "classname"            : prefix + basename,
                                  "smoothingScaleMethod" : smoothingScaleMethod})

#-------------------------------------------------------------------------------
# Provide a factory function to return the appropriate SVPH hydro.
#-------------------------------------------------------------------------------
def SVPH(dataBase,
         W,
         Q = None,
         cfl = 0.25,
         useVelocityMagnitudeForDt = False,
         useNewAccelerationMagnitudeForDt = False,
         compatibleEnergyEvolution = True,
         XSVPH = True,
         linearConsistent = False,
         generateVoid = True,
         densityUpdate = RigorousSumDensity,
         HUpdate = IdealH,
         fcentroidal = 0.0,
         fcellPressure = 0.0,
         xmin = (-1e100, -1e100, -1e100),
         xmax = ( 1e100,  1e100,  1e100),
         ASPH = False):

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    ndim = dataBase.nDim
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print("SVPH Error: you have provided both solid and fluid NodeLists, which is currently not supported.")
        print("            If you want some fluids active, provide SolidNodeList without a strength option specfied,")
        print("            which will result in fluid behaviour for those nodes.")
        raise RuntimeError("Cannot mix solid and fluid NodeLists.")

    # Decide on the hydro object.
    if ASPH:
        Constructor = eval("ASVPHHydro%id" % ndim)
    else:
        Constructor = eval("SVPHHydro%id" % ndim)

    # Artificial viscosity.
    if not Q:
        Cl = 1.0*(W.kernelExtent/2.0)
        Cq = 1.0*(W.kernelExtent/2.0)**2
        Q = eval("MonaghanGingoldViscosity%id(Clinear=%g, Cquadratic=%g)" % (ndim, Cl, Cq))

    # Build and return the thing.
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    result = Constructor(Q = Q,
                         W = W,
                         cfl = cfl,
                         useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                         useNewAccelerationMagnitudeForDt = useNewAccelerationMagnitudeForDt,
                         compatibleEnergyEvolution = compatibleEnergyEvolution,
                         XSVPH = XSVPH,
                         linearConsistent = linearConsistent,
                         generateVoid = generateVoid,
                         densityUpdate = densityUpdate,
                         HUpdate = HUpdate,
                         fcentroidal = fcentroidal,
                         fcellPressure = fcellPressure,
                         xmin = eval("Vector%id(%g, %g, %g)" % xmin),
                         xmax = eval("Vector%id(%g, %g, %g)" % xmax))
    result.Q = Q
    return result

#-------------------------------------------------------------------------------
# ASVPH factory.
#-------------------------------------------------------------------------------
def ASVPH(dataBase,
          W,
          Q = None,
          cfl = 0.25,
          useVelocityMagnitudeForDt = False,
          useNewAccelerationMagnitudeForDt = False,
          compatibleEnergyEvolution = True,
          XSVPH = True,
          linearConsistent = False,
          generateVoid = True,
          densityUpdate = RigorousSumDensity,
          HUpdate = IdealH,
          fcentroidal = 0.0,
          fcellPressure = 0.0,
          xmin = (-1e100, -1e100, -1e100),
          xmax = ( 1e100,  1e100,  1e100)):
    return SVPH(dataBase,
                W,
                Q = Q,
                cfl = cfl,
                useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                useNewAccelerationMagnitudeForDt = useNewAccelerationMagnitudeForDt,
                compatibleEnergyEvolution = compatibleEnergyEvolution,
                XSVPH = XSVPH,
                linearConsistent = linearConsistent,
                generateVoid = generateVoid,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                fcentroidal = fcentroidal,
                fcellPressure = fcellPressure,
                xmin = xmin,
                xmax = xmax,
                ASPH = True)
