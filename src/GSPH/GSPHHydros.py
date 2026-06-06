from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# GSPH convience wrapper function
#-------------------------------------------------------------------------------
def GSPH(dataBase,
         W,
         riemannSolver=None,
         specificThermalEnergyDiffusionCoefficient = 0.0,
         cfl = 0.25,
         gradientType = HydroAccelerationGradient,
         densityUpdate = IntegrateDensity,
         useVelocityMagnitudeForDt = False,
         compatibleEnergyEvolution = True,
         evolveTotalEnergy = False,
         XSPH = False,
         correctVelocityGradient = False,
         HUpdate = IdealH,
         epsTensile = 0.0,
         nTensile = 4.0,
         damageRelieveRubble = False,
         negativePressureInDamage = False,
         strengthInDamage = False,
         xmin = (-1e100, -1e100, -1e100),
         xmax = ( 1e100,  1e100,  1e100),
         ASPH = False,
         smoothingScaleMethod = None):

    assert densityUpdate in (RigorousSumDensity,IntegrateDensity)

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    ndim = dataBase.nDim
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print("GSPH Error: you have provided both solid and fluid NodeLists, which is currently not supported.")
        print("            If you want some fluids active, provide SolidNodeList without a strength option specfied,")
        print("            which will result in fluid behaviour for those nodes.")
        raise RuntimeError("Cannot mix solid and fluid NodeLists.")

    Constructor = eval("GSPH%id" % ndim)

    if riemannSolver is None:
        waveSpeedMethod = eval("DavisWaveSpeed%id()" % (ndim))
        slopeLimiter = eval("VanLeerLimiter%id()" % (ndim))
        linearReconstruction=True
        riemannSolver = eval("HLLC%id(slopeLimiter,waveSpeedMethod,linearReconstruction)" % (ndim))
   
    # Build the constructor arguments
    Vector = eval(f"Vector{ndim}d")
    kwargs = {"dataBase" : dataBase,
              "riemannSolver" : riemannSolver,
              "W" : W,
              "epsDiffusionCoeff" : specificThermalEnergyDiffusionCoefficient,
              "cfl" : cfl,
              "useVelocityMagnitudeForDt" : useVelocityMagnitudeForDt,
              "compatibleEnergyEvolution" : compatibleEnergyEvolution,
              "evolveTotalEnergy" : evolveTotalEnergy,
              "XSPH" : XSPH,
              "correctVelocityGradient" : correctVelocityGradient,
              "gradType" : gradientType,
              "densityUpdate" : densityUpdate,
              "epsTensile" : epsTensile,
              "nTensile" : nTensile,
              "xmin" : eval("Vector(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector(%g, %g, %g)" % xmax)}

    # Build and return the thing.
    result = Constructor(**kwargs)

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
# MFM convience wrapper function
#-------------------------------------------------------------------------------
def MFM(dataBase,
        W,
        riemannSolver=None,
        specificThermalEnergyDiffusionCoefficient = 0.0,
        cfl = 0.25,
        gradientType = HydroAccelerationGradient,
        densityUpdate = IntegrateDensity,
        useVelocityMagnitudeForDt = False,
        compatibleEnergyEvolution = True,
        evolveTotalEnergy = False,
        XSPH = False,
        correctVelocityGradient = False,
        HUpdate = IdealH,
        epsTensile = 0.0,
        nTensile = 4.0,
        damageRelieveRubble = False,
        negativePressureInDamage = False,
        strengthInDamage = False,
        xmin = (-1e100, -1e100, -1e100),
        xmax = ( 1e100,  1e100,  1e100),
        ASPH = False,
        smoothingScaleMethod = None):

    # for now we'll just piggy back off this enum
    assert densityUpdate in (RigorousSumDensity,IntegrateDensity)

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    ndim = dataBase.nDim
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print("MFM  Error: you have provided both solid and fluid NodeLists, which is currently not supported.")
        print("            If you want some fluids active, provide SolidNodeList without a strength option specfied,")
        print("            which will result in fluid behaviour for those nodes.")
        raise RuntimeError("Cannot mix solid and fluid NodeLists.")

    Constructor = eval("MFM%id" % ndim)

    if riemannSolver is None:
        waveSpeedMethod = eval("DavisWaveSpeed%id()" % (ndim))
        slopeLimiter = eval("VanLeerLimiter%id()" % (ndim))
        linearReconstruction=True
        riemannSolver = eval("HLLC%id(slopeLimiter,waveSpeedMethod,linearReconstruction)" % (ndim))
   
    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax

    kwargs = {"dataBase" : dataBase,
              "riemannSolver" : riemannSolver,
              "W" : W,
              "epsDiffusionCoeff" : specificThermalEnergyDiffusionCoefficient,
              "dataBase" : dataBase,
              "cfl" : cfl,
              "useVelocityMagnitudeForDt" : useVelocityMagnitudeForDt,
              "compatibleEnergyEvolution" : compatibleEnergyEvolution,
              "evolveTotalEnergy" : evolveTotalEnergy,
              "XSPH" : XSPH,
              "correctVelocityGradient" : correctVelocityGradient,
              "gradType" : gradientType,
              "densityUpdate" : densityUpdate,
              "epsTensile" : epsTensile,
              "nTensile" : nTensile,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}

    # Build and return the thing.
    result = Constructor(**kwargs)

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
# MFM convience wrapper function
#-------------------------------------------------------------------------------
def MFV(dataBase,
        W,
        riemannSolver=None,
        specificThermalEnergyDiffusionCoefficient = 0.0,
        cfl = 0.25,
        nodeMotionCoefficient = 0.2,
        nodeMotionType = NodeMotionType.Lagrangian,
        gradientType = HydroAccelerationGradient,
        densityUpdate = IntegrateDensity,
        useVelocityMagnitudeForDt = False,
        compatibleEnergyEvolution = True,
        evolveTotalEnergy = False,
        XSPH = False,
        correctVelocityGradient = False,
        HUpdate = IdealH,
        epsTensile = 0.0,
        nTensile = 4.0,
        damageRelieveRubble = False,
        negativePressureInDamage = False,
        strengthInDamage = False,
        xmin = (-1e100, -1e100, -1e100),
        xmax = ( 1e100,  1e100,  1e100),
        ASPH = False,
        RZ = False,
        smoothingScaleMethod = None):

    # for now we'll just piggy back off this enum
    assert densityUpdate in (RigorousSumDensity,IntegrateDensity)

    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    ndim = dataBase.nDim
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print("MFV  Error: you have provided both solid and fluid NodeLists, which is currently not supported.")
        print("            If you want some fluids active, provide SolidNodeList without a strength option specfied,")
        print("            which will result in fluid behaviour for those nodes.")
        raise RuntimeError("Cannot mix solid and fluid NodeLists.")

    Constructor = eval("MFV%id" % ndim)

    if riemannSolver is None:
        waveSpeedMethod = eval("DavisWaveSpeed%id()" % (ndim))
        slopeLimiter = eval("VanLeerLimiter%id()" % (ndim))
        linearReconstruction=True
        riemannSolver = eval("HLLC%id(slopeLimiter,waveSpeedMethod,linearReconstruction)" % (ndim))
   
    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax

    kwargs = {"dataBase" : dataBase,
              "riemannSolver" : riemannSolver,
              "W" : W,
              "epsDiffusionCoeff" : specificThermalEnergyDiffusionCoefficient,
              "dataBase" : dataBase,
              "cfl" : cfl,
              "useVelocityMagnitudeForDt" : useVelocityMagnitudeForDt,
              "compatibleEnergyEvolution" : compatibleEnergyEvolution,
              "evolveTotalEnergy" : evolveTotalEnergy,
              "XSPH" : XSPH,
              "correctVelocityGradient" : correctVelocityGradient,
              "nodeMotionCoefficient" : nodeMotionCoefficient,
              "nodeMotionType" : nodeMotionType,
              "gradType" : gradientType,
              "densityUpdate" : densityUpdate,
              "epsTensile" : epsTensile,
              "nTensile" : nTensile,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}

    # Build and return the thing.
    result = Constructor(**kwargs)

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

    #print(result.nodeMotionType)
    return result

#-------------------------------------------------------------------------------
# Provide shorthand names for ASPH
#-------------------------------------------------------------------------------
def AGSPH(*args, **kwargs):
    kwargs.update({"ASPH" : True})
    return GSPH(*args, **kwargs)

def AMFM(*args, **kwargs):
    kwargs.update({"ASPH" : True})
    return MFM(*args, **kwargs)

def AMFV(*args, **kwargs):
    kwargs.update({"ASPH" : True})
    return MFV(*args, **kwargs)
