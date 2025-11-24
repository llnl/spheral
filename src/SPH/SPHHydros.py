from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic SPHHydro pattern.
#-------------------------------------------------------------------------------
def SPH(W,
        WPi = None,
        WGrad = None,
        dataBase = None,
        Q = None,
        filter = None,
        cfl = 0.25,
        useVelocityMagnitudeForDt = False,
        compatibleEnergyEvolution = True,
        evolveTotalEnergy = False,
        gradhCorrection = True,
        XSPH = True,
        correctVelocityGradient = True,
        sumMassDensityOverAllNodeLists = True,
        densityUpdate = RigorousSumDensity,
        HUpdate = IdealH,
        epsTensile = 0.0,
        nTensile = 4.0,
        damageRelieveRubble = False,
        strengthInDamage = False,
        xmin = (-1e100, -1e100, -1e100),
        xmax = ( 1e100,  1e100,  1e100),
        etaMinAxis = 0.1,
        ASPH = False,
        smoothingScaleMethod = None):

    # Check if we're running solid or fluid hydro
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print("SPH Error: you have provided both solid and fluid NodeLists, which is currently not supported.")
        print("           If you want some fluids active, provide SolidNodeList without a strength option specfied,")
        print("           which will result in fluid behaviour for those nodes.")
        raise RuntimeError("Cannot mix solid and fluid NodeLists.")

    # Check for deprecated arguments
    if not filter is None:
        print("SPH DEPRECATION WARNING: filter is no longer used -- ignoring")

    # Pick the appropriate C++ constructor from dimensionality and coordinates
    ndim = dataBase.nDim
    if GeometryRegistrar.coords() == CoordinateType.Spherical:
        assert ndim == 1
        if nsolid > 0:
            constructor = SolidSphericalSPH
        else:
            constructor = SphericalSPH

        # Build the spherical kernels
        print("Constructing Spherical kernels...")
        W3S1 = SphericalKernel(W)
        W = W3S1
        if WPi:
            WPi3S1 = SphericalKernel(WPi)
            WPi = WPi3S1
        if WGrad:
            WGrad3S1 = SphericalKernel(WGrad)
            WGrad = WGrad3S1

    elif GeometryRegistrar.coords() == CoordinateType.RZ:
        assert ndim == 2
        if nsolid > 0:
            constructor = SolidSPHRZ
        else:
            constructor = SPHRZ

    else:
        if nsolid > 0:
            constructor = eval("SolidSPH%id" % ndim)
        else:
            constructor = eval("SPH%id" % ndim)

    # Fill out the set of kernels
    if WPi is None:
        WPi = W
    if WGrad is None:
        WGrad = W

    # Artificial viscosity.
    if not Q:
        Cl = 2.0*(dataBase.maxKernelExtent/2.0)
        Cq = 2.0*(dataBase.maxKernelExtent/2.0)**2
        if GeometryRegistrar.coords() == CoordinateType.Spherical:
            Q = eval("LimitedMonaghanGingoldViscosity%id(Clinear=%g, Cquadratic=%g, kernel=WPi.baseKernel1d)" % (ndim, Cl, Cq))
        else:
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
              "gradhCorrection" : gradhCorrection,
              "XSPH" : XSPH,
              "correctVelocityGradient" : correctVelocityGradient,
              "sumMassDensityOverAllNodeLists" : sumMassDensityOverAllNodeLists,
              "densityUpdate" : densityUpdate,
              "epsTensile" : epsTensile,
              "nTensile" : nTensile,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}

    if nsolid > 0:
        kwargs.update({"WGrad"                    : WGrad,
                       "damageRelieveRubble"      : damageRelieveRubble,
                       "strengthInDamage"         : strengthInDamage})

    # Build the SPH hydro
    result = constructor(**kwargs)

    # Add the Q as a sub-package (to run before the hydro)
    result.Q = Q
    result.prependSubPackage(Q)

    # Smoothing scale update
    if smoothingScaleMethod is None:
        WH = W.baseKernel1d if GeometryRegistrar.coords() == CoordinateType.Spherical else W
        if ASPH:
            if isinstance(ASPH, str) and ASPH.upper() == "CLASSIC":
                smoothingScaleMethod = eval(f"ASPHClassicSmoothingScale{ndim}d({HUpdate}, WH)")
            else:
                smoothingScaleMethod = eval(f"ASPHSmoothingScale{ndim}d({HUpdate}, WH)")
        else:
            smoothingScaleMethod = eval(f"SPHSmoothingScale{ndim}d({HUpdate}, WH)")
    result._smoothingScaleMethod = smoothingScaleMethod
    result.appendSubPackage(smoothingScaleMethod)

    # In spherical coordatinates, preserve our locally constructed spherical kernels
    if GeometryRegistrar.coords() == CoordinateType.Spherical:
        result._W3S1 = W
        result._WPi3S1 = WPi
        result._WGrad3S1 = WGrad

    # If we're using area-weighted RZ, store etaMinAxis
    if GeometryRegistrar.coords() == CoordinateType.RZ:
        result.etaMinAxis = etaMinAxis

    return result

#-------------------------------------------------------------------------------
# Provide shorthand names for ASPH
#-------------------------------------------------------------------------------
def ASPH(*args, **kwargs):
    kwargs.update({"ASPH" : True})
    return SPH(*args, **kwargs)
