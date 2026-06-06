from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic PSPHHydro pattern.
#-------------------------------------------------------------------------------
def PSPH(dataBase,
         W,
         WPi = None,
         Q = None,
         filter = None,
         cfl = 0.25,
         useVelocityMagnitudeForDt = False,
         compatibleEnergyEvolution = True,
         evolveTotalEnergy = False,
         XSPH = True,
         correctVelocityGradient = True,
         HopkinsConductivity = False,
         sumMassDensityOverAllNodeLists = True,
         densityUpdate = RigorousSumDensity,
         HUpdate = IdealH,
         xmin = (-1e100, -1e100, -1e100),
         xmax = ( 1e100,  1e100,  1e100),
         ASPH = False,
         smoothingScaleMethod = None):

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print("PSPH Warning: you have provided solid NodeLists, but PSPH currently does not have a solid option.")
        print("              The fluid limit will be provided for now.")

    # Check for deprecated arguments
    if not filter is None:
        print("PSPH DEPRECATION WARNING: filter is no longer used -- ignoring")

    # Pick the appropriate C++ constructor from dimensionality and coordinates
    ndim = dataBase.nDim
    constructor = eval("PSPH%id" % ndim)

    # Fill out the set of kernels
    if WPi is None:
        WPi = W

    # Artificial viscosity.
    if not Q:
        Cl = 2.0*(dataBase.maxKernelExtent/2.0)
        Cq = 2.0*(dataBase.maxKernelExtent/2.0)**2
        Q = eval("LimitedMonaghanGingoldViscosity%id(Clinear=%g, Cquadratic=%g, kernel=WPi)" % (ndim, Cl, Cq))

    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    kwargs = {"dataBase" : dataBase,
              "Q" : Q,
              "W" : W,
              "WPi" : WPi,
              "cfl" : cfl,
              "useVelocityMagnitudeForDt" : useVelocityMagnitudeForDt,
              "compatibleEnergyEvolution" : compatibleEnergyEvolution,
              "evolveTotalEnergy" : evolveTotalEnergy,
              "XSPH" : XSPH,
              "correctVelocityGradient" : correctVelocityGradient,
              "HopkinsConductivity" : HopkinsConductivity,
              "sumMassDensityOverAllNodeLists" : sumMassDensityOverAllNodeLists,
              "densityUpdate" : densityUpdate,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}

    # Build the thing
    result = constructor(**kwargs)
    result.Q = Q
    
    # Add the Q as a sub-package (to run before the hydro)
    result.prependSubPackage(Q)

    # Smoothing scale update
    if smoothingScaleMethod is None:
        if ASPH:
            if isinstance(ASPH, str) and ASPH.upper() == "CLASSIC":
                smoothingScaleMethod = eval(f"ASPHClassicSmoothingScale{ndim}d(HUpdate, W)")
            else:
                smoothingScaleMethod = eval(f"ASPHSmoothingScale{ndim}d(HUpdate, W)")
        else:
            smoothingScaleMethod = eval(f"SPHSmoothingScale{ndim}d(HUpdate, W)")
    result._smoothingScaleMethod = smoothingScaleMethod
    result.appendSubPackage(smoothingScaleMethod)

    return result

#-------------------------------------------------------------------------------
# Provide a shorthand APSPH factory.
#-------------------------------------------------------------------------------
def PASPH(*args, **kwargs):
    kwargs.update({"ASPH" : True})
    return PSPH(*args, **kwargs)
