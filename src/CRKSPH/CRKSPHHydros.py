from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic CRKSPHHydro pattern.
#-------------------------------------------------------------------------------
def CRKSPH(dataBase,
           W,
           Q = None,
           order = RKOrder.LinearOrder,
           filter = 0.0,
           cfl = 0.25,
           useVelocityMagnitudeForDt = False,
           compatibleEnergyEvolution = True,
           evolveTotalEnergy = False,
           XSPH = True,
           densityUpdate = RigorousSumDensity,
           HUpdate = IdealH,
           epsTensile = 0.0,
           nTensile = 4.0,
           damageRelieveRubble = False,
           ASPH = False,
           etaMinAxis = 0.1,
           crktype = "default",
           smoothingScaleMethod = None):

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print("CRKSPH Error: you have provided both solid and fluid NodeLists, which is currently not supported.")
        print("             If you want some fluids active, provide SolidNodeList without a strength option specfied,")
        print("             which will result in fluid behaviour for those nodes.")
        raise RuntimeError("Cannot mix solid and fluid NodeLists.")

    # Check for deprecated arguments
    if not filter is None:
        print("CRKSPH DEPRECATION WARNING: filter is no longer used -- ignoring")

    # Pick the appropriate C++ constructor from dimensionality and coordinates
    ndim = dataBase.nDim
    if GeometryRegistrar.coords() == CoordinateType.RZ:
        # RZ ----------------------------------------
        assert ndim == 2
        if nsolid > 0:
            constructor = SolidCRKSPHRZ
        else:
            constructor = CRKSPHRZ
    else:
        # Cartesian ---------------------------------
        crktype = crktype.lower()
        assert crktype in ("default", "variant")
        if nsolid > 0:
            constructor = eval("SolidCRKSPH%id" % ndim)
        else:
            if crktype == "variant":
                constructor = eval("CRKSPHVariant%id" % ndim)
            else:
                constructor = eval("CRKSPH%id" % ndim)

    # Artificial viscosity.
    if not Q:
        Cl = 2.0*(dataBase.maxKernelExtent/4.0)
        Cq = 1.0*(dataBase.maxKernelExtent/4.0)**2
        Q = eval("LimitedMonaghanGingoldViscosity%id(Clinear=%g, Cquadratic=%g, kernel=W)" % (ndim, Cl, Cq))

    # Build the constructor arguments
    kwargs = {"dataBase" : dataBase,
              "Q" : Q,
              "order" : order,
              "cfl" : cfl,
              "useVelocityMagnitudeForDt" : useVelocityMagnitudeForDt,
              "compatibleEnergyEvolution" : compatibleEnergyEvolution,
              "evolveTotalEnergy" : evolveTotalEnergy,
              "XSPH" : XSPH,
              "densityUpdate" : densityUpdate,
              "epsTensile" : epsTensile,
              "nTensile" : nTensile}

    if nsolid > 0:
        kwargs.update({"damageRelieveRubble" : damageRelieveRubble})

    # Build the thing.
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

    # If we're using area-weighted RZ, we need to reflect from the axis
    if GeometryRegistrar.coords() == CoordinateType.RZ:
        result.zaxisBC = AxisBoundaryRZ(etaMinAxis)
        result.appendBoundary(result.zaxisBC)

    return result

#-------------------------------------------------------------------------------
# ACRKSPH
#-------------------------------------------------------------------------------
def ACRKSPH(*args, **kwargs):
    kwargs.update({"ASPH" : True})
    return CRKSPH(*args, **kwargs)
